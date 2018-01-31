"""
```
decompose_forecast_one_draw(model, df_new, df_old, params_new, params_old, k, cond_type)
```

Returns reestimation, state_updates, realized_shocks, data_revisions, overall_change
"""
function decompose_forecast_one_draw(model::AbstractModel, df_new::DataFrame, df_old::DataFrame,
                                     params_new::Vector{Float64}, params_old::Vector{Float64}, k::Int64;
                                     cond_type::Symbol = :none)

    # Update model
    DSGE.update!(model, params_new)
    sys_new = compute_system(model)
    DSGE.update!(model, params_old)
    sys_old = compute_system(model)

    # Filter and smooth
    kal_new = DSGE.filter(model, df_new, sys_new; cond_type = cond_type) # |T, new para
    kal_mix = DSGE.filter(model, df_old, sys_new; cond_type = cond_type) # |T-K, new para
    kal_old = DSGE.filter(model, df_old, sys_old; cond_type = cond_type) # |T-k, old para

    histstates, histshocks, _, _ = DSGE.smooth(model, df_new, sys_new, kal_new; draw_states = false, cond_type = cond_type)

    # Compute decomposition
    n_obs    = n_observables(model)
    n_shocks = n_shocks_exogenous(model)
    n_states = n_states_augmented(model)

    reestimation    = zeros(n_obs, forecast_horizons(model))
    data_revisions  = zeros(n_obs, forecast_horizons(model))
    state_updates   = zeros(n_obs, forecast_horizons(model))
    realized_shocks = zeros(n_obs, forecast_horizons(model), n_shocks)
    overall_change  = zeros(n_obs, forecast_horizons(model))

    Z, T, R, D, C = sys_old[:ZZ], sys_old[:TTT], sys_old[:RRR], sys_old[:DD], sys_old[:CCC]
    Z_new, T_new, R_new, D_new, C_new = sys_new[:ZZ], sys_new[:TTT], sys_new[:RRR], sys_new[:DD], sys_new[:CCC]

    c = cond_type in [:semi, :full] ? 1 : 0
    partial_sum(M::Matrix, k::Int) = k <= 0 ? eye(size(M)...) : M^k + partial_sum(M, k-1)

    state_component(h)    = Z_new * T_new^(h+k) * (histstates[:,end-k] - kal_new[:filt][:,end-k])
    revision_component(h) = Z_new * T_new^(h+k) * (kal_new[:filt][:,end-k] - kal_mix[:filt][:,end])
    para_component(h)     = Z_new * (T_new^(h+k) * kal_mix[:filt][:,end] + partial_sum(T_new,h+k) * C_new) - Z * (T^(h+k) * kal_old[:filt][:,end] + partial_sum(T,h+k) * C) + (D_new - D)
    all_components(h)     = Z_new * (T_new^h * kal_new[:zend] + C_new) - Z * (T^(h+k) * kal_old[:zend] + C) + (D_new - D)
    function shock_component(h::Int)
        total = zeros(n_obs, 1, n_shocks)
        for s in 1:n_shocks
            shock_total = zeros(n_states, 1)
            shock = zeros(n_shocks, 1)
            for i in 0:k-1
                shock[s,1] = histshocks[s,end-i]
                shock_total += T_new^(h+i) * R_new * shock
            end
            total[:,1,s] = Z_new * shock_total
        end
        return total
    end

    for t in 1:forecast_horizons(model)
        reestimation[:,t]      = para_component(t-c)
        state_updates[:,t]     = state_component(t-c)
        realized_shocks[:,t,:] = shock_component(t-c)
        data_revisions[:,t]    = revision_component(t-c)
        overall_change[:,t]    = all_components(t-c)
    end

    return reestimation, state_updates, realized_shocks, data_revisions, overall_change

end


"""
```
decompose_forecast_changes(m_new, m_old, input_type, cond_type;
                           verbose = :low, transform = false,
                           format_as_DataFrame= true)

decompose_forecast_changes(m_new, m_old, df_new, df_old, input_type, cond_type;
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
                               verbose = verbose, transform = transform,
                               format_as_DataFrame = format_as_DataFrame)
end

function decompose_forecast_changes(m_new::AbstractModel, m_old::AbstractModel,
                                    df_new::DataFrame, df_old::DataFrame,
                                    input_type::Symbol, cond_type::Symbol;
                                    verbose::Symbol = :low, transform::Bool = false,
                                    format_as_DataFrame = true)

    # calculate offset
    new_start = date_forecast_start(m_new)
    old_start = date_forecast_start(m_old)
    quarter_diff = DSGE.quarter_date_to_number(new_start) - DSGE.quarter_date_to_number(old_start)
    offset = Int(4 * quarter_diff)
    @assert offset >= 0 "The vintage of m_new must not be older than the vintage of m_old"

    # Initialize results matrices
    reestimation    = zeros(n_observables(m_new), forecast_horizons(m_new))
    state_updates   = zeros(n_observables(m_new), forecast_horizons(m_new))
    realized_shocks = zeros(n_observables(m_new), forecast_horizons(m_new), n_shocks_exogenous(m_new))
    data_revisions  = zeros(n_observables(m_new), forecast_horizons(m_new))
    overall_change  = zeros(n_observables(m_new), forecast_horizons(m_new))

    decompose(params) = decompose_forecast_one_draw(m_new, df_new, df_old, params[1], params[2], offset;
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
            reestimation    = reduce(+, reestimation, map(x->x[1], forecast_outputs))
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
        reestimation    = DSGE.quartertoannual(reestimation)
        state_updates   = DSGE.quartertoannual(state_updates)
        realized_shocks = DSGE.quartertoannual(realized_shocks)
        data_revisions  = DSGE.quartertoannual(data_revisions)
        overall_change  = DSGE.quartertoannual(overall_change)
    end

    if format_as_DataFrame
        n_periods  = forecast_horizons(m_new)
        start_date = date_forecast_start(m_new)
        end_date   = iterate_quarters(start_date, n_periods - 1)
        var_names = sort(m_new.observables, by=x->m_new.observables[x]).keys

        function forecast_mat_to_df(mat, names)
            df = DataFrame()
            df[:date] = DSGE.quarter_range(start_date, end_date)

            for (i,name) in enumerate(names)
                df[name] = vec(mat[i,:])
            end
            return df
        end

        reestimation    = forecast_mat_to_df(reestimation, var_names)
        state_updates   = forecast_mat_to_df(state_updates, var_names)
        data_revisions  = forecast_mat_to_df(data_revisions, var_names)
        overall_change  = forecast_mat_to_df(overall_change, var_names)

        if false
            realized_shocks = forecast_mat_to_df(sum(realized_shocks, 3), df, var_names)
        else
            realized_shocks_new = Vector{DataFrames.DataFrame}(n_shocks_exogenous(m_new))
            for s in 1:n_shocks_exogenous(m_new)
                new_var_names = map(x -> Symbol(x, :_, m_new.exogenous_shocks.keys[s]), var_names)
                realized_shocks_new[s] = forecast_mat_to_df(realized_shocks[:,:,s], new_var_names)
            end

            realized_shocks = realized_shocks_new
        end
    end

return reestimation, state_updates, realized_shocks, data_revisions, overall_change
end


function collapse_shocks(shocks::Vector{DataFrames.DataFrame}, m::AbstractModel; cond_type::Symbol = :none)
    n_shocks = n_shocks_exogenous(m)
    var_names = sort(m.observables, by=x->m.observables[x]).keys

    df = DataFrame()
    n_periods = size(shocks[1],1)
    start_date = date_forecast_start(m_new)
    end_date = iterate_quarters(start_date, n_periods - 1)
    df[:date] = DSGE.quarter_range(start_date, end_date)

    for (i,var) in enumerate(var_names)
        total = Vector{Float64}(size(shocks[1],1))
        for s in 1:n_shocks
            total += Vector(shocks[s][:,i+1])
        end
        df[var] = total
    end
    return df
end

function collapse_shocks(shocks::Array{Float64,3})
    return squeeze(sum(shocks,3),3)
end
