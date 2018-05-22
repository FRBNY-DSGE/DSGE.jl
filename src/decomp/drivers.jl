# TODO:
# - Maybe: handle h < 0 (new history vs. old history)
# - Maybe: break out checking into own functions
# - Plotting

function decompose_forecast(m_new::AbstractModel, m_old::AbstractModel,
                            df_new::DataFrame, df_old::DataFrame,
                            input_type::Symbol, cond_new::Symbol, cond_old::Symbol,
                            classes::Vector{Symbol};
                            hs = 1:forecast_horizons(m_old),
                            verbose::Symbol = :low, kwargs...)
    # Get output file names
    decomp_output_files = get_decomp_output_files(m_new, m_old, input_type, cond_new, cond_old, classes)

    if DSGE.VERBOSITY[verbose] >= DSGE.VERBOSITY[:low]
        info("Decomposing forecast...")
        println("Start time: $(now())")
    end

    # Single-draw forecasts
    if input_type in [:mode, :mean, :init]

        params_new = load_draws(m_new, input_type, verbose = verbose)
        params_old = load_draws(m_old, input_type, verbose = verbose)

        decomps = decompose_forecast(m_new, m_old, df_new, df_old, params_new, params_old,
                                     cond_new, cond_old, classes; hs = hs, kwargs...)
        write_forecast_decomposition(m_new, m_old, input_type, classes, hs, decomp_output_files, decomps,
                                     verbose = verbose)

    # Multiple-draw forecasts
    elseif input_type == :full

        block_inds, block_inds_thin = DSGE.forecast_block_inds(m_new, input_type)
        nblocks = length(block_inds)
        total_forecast_time = 0.0

        for block = 1:nblocks
            if DSGE.VERBOSITY[verbose] >= DSGE.VERBOSITY[:low]
                println()
                info("Decomposing block $block of $nblocks...")
            end
            tic()

            # Get to work!
            params_new = load_draws(m_new, input_type, block_inds[block], verbose = verbose)
            params_old = load_draws(m_old, input_type, block_inds[block], verbose = verbose)

            mapfcn = use_parallel_workers(m_new) ? pmap : map
            decomps = mapfcn((param_new, param_old) ->
                             decompose_forecast(m_new, m_old, df_new, df_old, param_new, param_old,
                                                cond_new, cond_old, classes; hs = hs, kwargs...),
                             params_new, params_old)

            # Assemble outputs from this block and write to file
            decomps = convert(Vector{Dict{Symbol, Array{Float64}}}, decomps)
            decomps = DSGE.assemble_block_outputs(decomps)
            write_forecast_decomposition(m_new, m_old, input_type, classes, hs, decomp_output_files, decomps,
                                         block_number = Nullable(block), block_inds = block_inds_thin[block],
                                         verbose = verbose)
            gc()

            # Calculate time to complete this block, average block time, and
            # expected time to completion
            if DSGE.VERBOSITY[verbose] >= DSGE.VERBOSITY[:low]
                block_time = toq()
                total_forecast_time += block_time
                total_forecast_time_min     = total_forecast_time/60
                expected_time_remaining     = (total_forecast_time/block)*(nblocks - block)
                expected_time_remaining_min = expected_time_remaining/60

                println("\nCompleted $block of $nblocks blocks.")
                println("Total time elapsed: $total_forecast_time_min minutes")
                println("Expected time remaining: $expected_time_remaining_min minutes")
            end
        end # of loop through blocks

    else
        error("Invalid input_type: $input_type. Must be in [:mode, :mean, :init, :full]")
    end

    if DSGE.VERBOSITY[verbose] >= DSGE.VERBOSITY[:low]
        println("\nForecast decomposition complete: $(now())")
    end
end

function decompose_forecast(m_new::AbstractModel, m_old::AbstractModel,
                            df_new::DataFrame, df_old::DataFrame,
                            params_new::Vector{Float64}, params_old::Vector{Float64},
                            cond_new::Symbol, cond_old::Symbol, classes::Vector{Symbol};
                            hs = 1:forecast_horizons(m_old),
                            check::Bool = false, atol::Float64 = 1e-8)

    # Compute numbers of periods
    # h_cond and k_cond are h and k adjusted for (differences in) conditioning
    T0, T, k, T1_new, T1_old, k_cond =
        decomposition_periods(m_new, m_old, cond_new, cond_old)

    if check
        @assert size(df_new, 1) == T0 + T + T1_new
        @assert size(df_old, 1) == T0 + T + T1_old - k
    end

    # Update parameters, filter, and smooth
    sys_new, sys_old, kal_new_new, kal_new_old, kal_old_old, s_tgT, ϵ_tgT =
        prepare_decomposition!(m_new, m_old, df_new, df_old, params_new, params_old,
                               cond_new, cond_old)

    # Initialize output dictionary
    decomp = Dict{Symbol, Array{Float64}}()
    for class in classes
        ZZ, _, _, _, _ = class_system_matrices(sys_new, class)
        nobs = size(ZZ, 1)
        for comp in [:state, :shock, :data, :param, :total]
            output_var = Symbol(:decomp, comp, class)
            decomp[output_var] = zeros(size(ZZ, 1), length(hs))
        end
    end

    # Compute whole sequence of historical and forecasted observables if desired
    if check
        y_new = Dict{Symbol, Matrix{Float64}}()
        y_old = Dict{Symbol, Matrix{Float64}}()
        for class in classes
            # y^{new,new}_{t|T}, t = 1:T+H
            y_new[class] = compute_history_and_forecast(m_new, df_new, class, cond_type = cond_new)

            # y^{old,old}_{t|T-k}, t = 1:T-k+H
            y_old[class] = compute_history_and_forecast(m_old, df_old, class, cond_type = cond_old)
        end
    end

    for (i, h) in enumerate(hs)
        # Adjust h for (differences in) conditioning
        h_cond = h - T1_old
        @assert h_cond >= 0

        for class in classes
            # Decompose into components
            decomp[Symbol(:decompstate, class)][:, i], decomp[Symbol(:decompshock, class)][:, i] =
                decompose_states_shocks(sys_new, kal_new_new, s_tgT, ϵ_tgT,
                                        class, T0, k_cond, h_cond, check = check)
            decomp[Symbol(:decompdata, class)][:, i] =
                decompose_data_revisions(sys_new, kal_new_new, kal_new_old,
                                         class, T0, k_cond, h_cond, check = check)
            decomp[Symbol(:decompparam, class)][:, i] =
                decompose_param_reest(sys_new, sys_old, kal_new_old, kal_old_old,
                                      class, T0, k_cond, h_cond, check = check)

            # Compute total difference
            decomp[Symbol(:decomptotal, class)][:, i] =
                decomp[Symbol(:decompstate, class)][:, i] + decomp[Symbol(:decompshock, class)][:, i] +
                decomp[Symbol(:decompdata, class)][:, i]  + decomp[Symbol(:decompparam, class)][:, i]

            if check
                # Total difference = y^{new,new}_{T-k+h|T} - y^{old,old}_{T-k+h|T-k}
                y_new_Tmkph_T   = y_new[class][:, T-k+h]
                y_old_Tmkph_Tmk = y_old[class][:, T-k+h]
                exp_total = y_new_Tmkph_T - y_old_Tmkph_Tmk
                try
                    @assert isapprox(exp_total, decomp[Symbol(:decomptotal, class)][:, i], atol = atol)
                catch ex
                    @show decomp[Symbol(:decomptotal, class)][:, i]
                    @show exp_total
                    throw(ex)
                end
            end
        end
    end

    return decomp
end

function decomposition_periods(m_new::AbstractModel, m_old::AbstractModel,
                               cond_new::Symbol, cond_old::Symbol)
    # Number of presample periods T0 must be the same
    T0 = n_presample_periods(m_new)
    @assert n_presample_periods(m_old) == T0

    # New model has T main-sample periods
    # Old model has T-k main-sample periods
    T = n_mainsample_periods(m_new)
    k = DSGE.subtract_quarters(date_forecast_start(m_new), date_forecast_start(m_old))
    @assert k >= 0

    # Number of conditional periods T1 may differ
    T1_new = cond_new == :none ? 0 : n_conditional_periods(m_new)
    T1_old = cond_old == :none ? 0 : n_conditional_periods(m_old)

    # Adjust k for (differences in) conditioning
    k_cond = k + T1_new - T1_old

    return T0, T, k, T1_new, T1_old, k_cond
end

function prepare_decomposition!(m_new::AbstractModel, m_old::AbstractModel,
                                df_new::DataFrame, df_old::DataFrame,
                                params_new::Vector{Float64}, params_old::Vector{Float64},
                                cond_new::Symbol, cond_old::Symbol)
    # Check models well-formed
    @assert typeof(m_new) == typeof(m_old)

    # Update parameters
    DSGE.update!(m_new, params_new)
    DSGE.update!(m_old, params_old)
    sys_new = compute_system(m_new)
    sys_old = compute_system(m_old)

    # Filter and smooth
    kal_new_new = DSGE.filter(m_new, df_new, sys_new, cond_type = cond_new)
    kal_new_old = DSGE.filter(m_new, df_old, sys_new, cond_type = cond_old)
    kal_old_old = DSGE.filter(m_old, df_old, sys_old, cond_type = cond_old)
    s_tgT, ϵ_tgT = smooth(m_new, df_new, sys_new, cond_type = cond_new, draw_states = false)

    return sys_new, sys_old, kal_new_new, kal_new_old, kal_old_old, s_tgT, ϵ_tgT
end

function decompose_states_shocks(sys_new::System, kal_new::DSGE.Kalman,
                                 s_tgT::Matrix{Float64}, ϵ_tgT::Matrix{Float64},
                                 class::Symbol, T0::Int, k::Int, h::Int;
                                 check::Bool = false, atol::Float64 = 1e-8)
    # New parameters, new data
    # State and shock components = y_{T-k+h|T} - y_{T-k+h|T-k}
    #     = Z [ T^h (s_{T-k|T} - s_{T-k|T-k}) + sum_{j=1}^(min(k,h)) ( T^(h-j) R ϵ_{T-k+j|T} ) ]
    #     = state component + shock component
    ZZ, DD, TTT, RRR, CCC = class_system_matrices(sys_new, class)
    s_tgt = kal_new[:s_filt][:, (T0+1):end] # s_{t|t}

    # State component = Z T^h (s_{T-k|T} - s_{T-k|T-k})
    s_Tmk_T   = s_tgT[:, end-k] # s_{T-k|T}
    s_Tmk_Tmk = s_tgt[:, end-k] # s_{T-k|T-k}
    state_comp = ZZ * TTT^h * (s_Tmk_T - s_Tmk_Tmk)

    # Shock component = Z sum_{j=1}^(min(k,h)) ( T^(h-j) R ϵ_{T-k+j|T} )
    shock_sum = zeros(size(s_Tmk_T))
    for j = 1:min(k, h)
        shock_sum += TTT^(h-j) * RRR * ϵ_tgT[:, end-k+j]
    end
    shock_comp = ZZ * shock_sum

    # Return
    if check
        @assert size(s_tgt, 2) == size(s_tgT, 2) == size(ϵ_tgT, 2)
        @assert isapprox(s_tgt[:, end], s_tgT[:, end], atol = atol)

        y_Tmkph_T = if h <= k
            ZZ * (s_tgT[:, end-k+h] + CCC) + DD # history
        else
            ZZ * (TTT^(h-k) * s_tgt[:, end] + CCC) + DD # forecast
        end
        y_Tmkph_Tmk = ZZ * (TTT^h * s_tgt[:, end-k] + CCC) + DD
        state_shock_comp = y_Tmkph_T - y_Tmkph_Tmk

        try
            @assert isapprox(state_shock_comp, state_comp + shock_comp, atol = atol)
        catch ex
            @show state_comp + shock_comp
            @show state_shock_comp
            throw(ex)
        end
    end
    return state_comp, shock_comp
end

function decompose_data_revisions(sys_new::System, kal_new::DSGE.Kalman, kal_old::DSGE.Kalman,
                                  class::Symbol, T0::Int, k::Int, h::Int;
                                  check::Bool = false, atol::Float64 = 1e-8)
    # New parameters, new and old data
    # Data revision component = y^{new}_{T-k+h|T-k} - y^{old}_{T-k+h|T-k}
    #     = Z T^h ( s^{new}_{T-k|T-k} - s^{old}_{T-k|T-k} )
    ZZ, _, TTT, _, _ = class_system_matrices(sys_new, class)

    s_new_tgt = kal_new[:s_filt][:, (T0+1):end] # s^{new}_{t|t}, t = 1:T
    s_old_tgt = kal_old[:s_filt][:, (T0+1):end] # s^{old}_{t|t}, t = 1:T-k
    try
        check && @assert size(s_new_tgt, 2) == size(s_old_tgt, 2) + k
    catch ex
        @show size(s_new_tgt, 2)
        @show size(s_old_tgt, 2)
        @show k
        throw(ex)
    end

    # Return
    data_comp = ZZ * TTT^h * (s_new_tgt[:, end-k] - s_old_tgt[:, end])
    return data_comp
end

function decompose_param_reest(sys_new::System, sys_old::System,
                               kal_new::DSGE.Kalman, kal_old::DSGE.Kalman,
                               class::Symbol, T0::Int, k::Int, h::Int;
                               check::Bool = false, atol::Float64 = 1e-8)
    # New and old parameters, old data
    # Parameter re-estimation component = y^{new}_{T-k+h|T-k} - y^{old}_{T-k+h|T-k}
    ZZ_new, DD_new, TTT_new, _, CCC_new = class_system_matrices(sys_new, class)
    ZZ_old, DD_old, TTT_old, _, CCC_old = class_system_matrices(sys_old, class)

    s_new_tgt = kal_new[:s_filt][:, (T0+1):end] # s^{new}_{t|t}, t = 1:T-k
    s_old_tgt = kal_old[:s_filt][:, (T0+1):end] # s^{old}_{t|t}, t = 1:T-k
    if check
        @assert size(s_new_tgt, 2) == size(s_old_tgt, 2)
    end

    # y^{new}_{T-k+h|T-k} = Z^{new} (T^{new}^h s^{new}_{T-k|T-k} + C^{new}) + D^{new}
    s_new_Tmk_Tmk = s_new_tgt[:, end]
    y_new_Tmkph_Tmk = ZZ_new * (TTT_new^h * s_new_Tmk_Tmk + CCC_new) + DD_new

    # y^{old}_{T-k+h|T-k} = Z^{old} (T^{old}^h s^{old}_{T-k|T-k} + C^{old}) + D^{old}
    s_old_Tmk_Tmk = s_old_tgt[:, end]
    y_old_Tmkph_Tmk = ZZ_old * (TTT_old^h * s_old_Tmk_Tmk + CCC_old) + DD_old

    # Return
    param_comp = y_new_Tmkph_Tmk - y_old_Tmkph_Tmk
    return param_comp
end