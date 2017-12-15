"""
```
decompose_forecast_one_draw(model, df_new, df_old, params)
```

Returns reestimation, state_updates, realized_shocks, data_revisions, overall_change
"""
function decompose_forecast_one_draw(model::AbstractModel, df_new::DataFrame, df_old::DataFrame,
                                     params_new::Vector{Float64}, params_old::Vector{Float64};
                                     cond_type::Symbol = :none)

    # Update model
    DSGE.update!(model, params_new)
    sys_new = compute_system(model)
    DSGE.update!(model, params_old)
    sys_old = compute_system(model)

    # Filter and smooth
    kal_new = DSGE.filter(model, df_new, sys_new; cond_type = cond_type) # |T, new para
    kal_mix = DSGE.filter(model, df_new, sys_old; cond_type = cond_type) # |T, old para
    kal_old = DSGE.filter(model, df_old, sys_old; cond_type = cond_type) # |T-1, old para

    histstates_mix, histshocks_mix, _, _ = DSGE.smooth(model, df_new, sys_old, kal_mix; draw_states = false)

    # Compute decomposition
    reestimation    = zeros(n_observables(model), forecast_horizons(model))
    data_revisions  = zeros(n_observables(model), forecast_horizons(model))
    state_updates   = zeros(n_observables(model), forecast_horizons(model))
    realized_shocks = zeros(n_observables(model), forecast_horizons(model))
    overall_change  = zeros(n_observables(model), forecast_horizons(model))

    Z, T, R = sys_old[:ZZ], sys_old[:TTT], sys_old[:RRR]
    Z_new, T_new = sys_new[:ZZ], sys_new[:TTT]

    para_component(j)     = Z_new * T_new^j * kal_new[:zend] - Z * T^j * kal_mix[:zend]
    state_component(j)    = Z * T^(j+1) * (histstates_mix[:,end-1] - kal_mix[:filt][:,end-1])
    shock_component(j)    = Z * T^j * R * histshocks_mix[:,end]
    revision_component(j) = Z * T^(j+1) * (kal_mix[:filt][:,end-1] - kal_old[:zend])
    all_components(j)     = Z_new * T_new^j * kal_new[:zend] - Z * T^(j+1) * kal_old[:zend]

    for t in 1:forecast_horizons(model)
        reestimation[:,t]    = para_component(t)
        state_updates[:,t]   = state_component(t)
        realized_shocks[:,t] = shock_component(t)
        data_revisions[:,t]  = revision_component(t)
        overall_change[:,t]  = all_components(t)
    end

    return reestimation, state_updates, realized_shocks, data_revisions, overall_change

end


"""
```
decompose_forecast_changes(m_new, m_old, df_new, df_old, input_type, cond_type;
                           verbose = :low, transform = false,
                           format_as_DataFrame= true)

decompose_forecast_changes(m_new, m_old, input_type, cond_type;
                           verbose = :low, transform = false,
                           format_as_DataFrame= true)
```

Decompose the difference between forecasts in `yr`:`quarter` and the subsequent quarter,
(e.g. forecasts in 2010Q1 and 2010Q2) into components attributable to updating the states,
new shocks, and data revisions. Returns changes due to state updates, realized shocks,
and data revisions
"""
function decompose_forecast_changes(m_new::AbstractModel, m_old::AbstractModel,
                                    input_type::Symbol, cond_type::Symbol;
                                    verbose::Symbol = :low, transform::Bool = false,
                                    format_as_DataFrame = true)

    # Load data
    df_new = load_data(m_new; cond_type = cond_type)
    df_old = load_data(m_old; cond_type = cond_type)

    decompose_forecast_changes(m_new, m_old, df_new, df_old, input_type, cond_type;
                               verbose = verbose, transform = transformd,
                               format_as_DataFrame = format_as_DataFrame)
end

function decompose_forecast_changes(m_new::AbstractModel, m_old::AbstractModel,
                                    df_new::DataFrame, df_old::DataFrame,
                                    input_type::Symbol, cond_type::Symbol;
                                    verbose::Symbol = :low, transform::Bool = false,
                                    format_as_DataFrame = true)

    # Initialize results matrices
    reestimation    = zeros(n_observables(m_new), forecast_horizons(m_new))
    state_updates   = zeros(n_observables(m_new), forecast_horizons(m_new))
    realized_shocks = zeros(n_observables(m_new), forecast_horizons(m_new))
    data_revisions  = zeros(n_observables(m_new), forecast_horizons(m_new))
    overall_change  = zeros(n_observables(m_new), forecast_horizons(m_new))

    decompose(params) = decompose_forecast_one_draw(m_new, df_new, df_old, params[1], params[2];
                                                    cond_type = cond_type)

    # Choose parameters depending on input_type and run decomposition
    if input_type == :init
        init_params = map(x->x.value, m_new.parameters)
        reestimation, state_updates, realized_shocks, data_revisions, overall_change = decompose((init_params, init_params))

    elseif input_type in [:mode, :mean]
        new_params = load_draws(m_new, input_type; verbose = verbose)
        old_params = load_draws(m_old, input_type; verbose = verbose)
        reestimation, state_updates, realized_shocks, data_revisions, overall_change = decompose((new_params, old_params))

    elseif input_type == :full

        # Block info
        block_inds, block_inds_thin = DSGE.forecast_block_inds(m_new, input_type)
        nblocks = length(block_inds)
        start_block = isnull(forecast_start_block(m_new)) ? 1 : get(forecast_start_block(m_new))

        # Info needed for printing progress
        total_forecast_time = 0.0

        for block = start_block:nblocks
            println()
            info("Forecasting block $block of $nblocks...")

            tic()

            # How we match old and new parameters does not matter since we average
            new_params = load_draws(m_new, input_type, block_inds[block]; verbose = verbose)
            old_params = load_draws(m_old, input_type, block_inds[block]; verbose = verbose)

            # Do the work
            forecast_outputs = pmap(decompose, zip(new_params, old_params))

            # If some element of forecast_outputs is a RemoteException, rethrow the exception
            ind_ex = findfirst(x -> isa(x, RemoteException), forecast_outputs)
            if ind_ex > 0
                ex = forecast_outputs[ind_ex].captured
                throw(ex)
            end

            # Add results from this block to results matrices
            reestimation    = reduce(+, state_updates, map(x->x[1], forecast_outputs))
            state_updates   = reduce(+, state_updates, map(x->x[2], forecast_outputs))
            realized_shocks = reduce(+, realized_shocks, map(x->x[3], forecast_outputs))
            data_revisions  = reduce(+, data_revisions, map(x->x[4], forecast_outputs))
            overall_change  = reduce(+, overall_change, map(x->x[5], forecast_outputs))

            gc()

            # Calculate time to complete this block, average block time, and
            # expected time to completion
            block_time = toq()
            total_forecast_time += block_time
            total_forecast_time_min     = total_forecast_time/60
            blocks_elapsed              = block - start_block + 1
            expected_time_remaining     = (total_forecast_time/blocks_elapsed)*(nblocks - block)
            expected_time_remaining_min = expected_time_remaining/60

            println("\nCompleted $block of $nblocks blocks.")
            println("Total time elapsed: $total_forecast_time_min minutes")
            println("Expected time remaining: $expected_time_remaining_min minutes")

        end # of loop through blocks

        # Convert sum to average
        n_draws = get_setting(m_new, :forecast_block_size) * nblocks
        reestimation    = reestimation / n_draws
        state_updates   = state_updates / n_draws
        realized_shocks = realized_shocks / n_draws
        data_revisions  = data_revisions / n_draws
        overall_change  = overall_change / n_draws

    else
        error("Invalid input_type")
    end

    if transform
        # quartertoannual works as a first approximation (ignoring per-capita
        # concerns) to most reverse transformations
        reestimation    = quartertoannual(reestimation)
        state_updates   = quartertoannual(state_updates)
        realized_shocks = quartertoannual(realized_shocks)
        data_revisions  = quartertoannual(data_revisions)
        overall_change  = quartertoannual(overall_change)
    end

    if format_as_DataFrame
        function forecast_mat_to_df(mat, init_df, var_names)
            df = deepcopy(init_df)
            for (i,var_name) in enumerate(var_names)
                df[var_name] = vec(mat[i,:])
            end
            return df
        end

        df = DataFrame()
        n_periods = size(state_updates, 2)
        start_date = date_forecast_start(m_new)
        end_date = iterate_quarters(start_date, n_periods - 1)
        df[:date] = DSGE.quarter_range(start_date, end_date)
        var_names = sort(m_new.observables, by=x->m_new.observables[x]).keys

        reestimation    = forecast_mat_to_df(reestimation, df, var_names)
        state_updates   = forecast_mat_to_df(state_updates, df, var_names)
        realized_shocks = forecast_mat_to_df(realized_shocks, df, var_names)
        data_revisions  = forecast_mat_to_df(data_revisions, df, var_names)
        overall_change  = forecast_mat_to_df(overall_change, df, var_names)

    end


return reestimation, state_updates, realized_shocks, data_revisions, overall_change
end
