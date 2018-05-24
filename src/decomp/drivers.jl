# TODO:
# - Maybe: handle h < 0 (new history vs. old history)
# - Plotting

function decompose_forecast(m_new::M, m_old::M, df_new::DataFrame, df_old::DataFrame,
                            input_type::Symbol, cond_new::Symbol, cond_old::Symbol,
                            classes::Vector{Symbol};
                            hs = 1:forecast_horizons(m_old),
                            verbose::Symbol = :low, kwargs...) where M<:AbstractModel
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

function decompose_forecast(m_new::M, m_old::M, df_new::DataFrame, df_old::DataFrame,
                            params_new::Vector{Float64}, params_old::Vector{Float64},
                            cond_new::Symbol, cond_old::Symbol, classes::Vector{Symbol};
                            hs = 1:forecast_horizons(m_old),
                            check::Bool = false, atol::Float64 = 1e-8) where M<:AbstractModel

    # Compute numbers of periods
    # k_cond is h adjusted for (differences in) conditioning
    T0, T, k, T1_new, T1_old, k_cond =
        decomposition_periods(m_new, m_old, cond_new, cond_old)

    @assert size(df_new, 1) == T0 + T + T1_new
    @assert size(df_old, 1) == T0 + T + T1_old - k

    # Update parameters, filter, and smooth
    sys_new, sys_old, s_new_new_tgt, s_new_new_tgT, ϵ_new_new_tgT, s_old_new_Tmk_Tmk, s_old_old_Tmk_Tmk =
        prepare_decomposition!(m_new, m_old, df_new, df_old, params_new, params_old,
                               cond_new, cond_old, k_cond; atol = atol)

    s_new_new_Tmk_Tmk = s_new_new_tgt[:, end-k_cond]
    s_new_new_Tmk_T   = s_new_new_tgT[:, end-k_cond]

    # Initialize output dictionary
    decomp = Dict{Symbol, Array{Float64}}()
    for class in classes
        ZZ, _ = class_measurement_matrices(sys_new, class)
        Ny    = size(ZZ, 1)
        for comp in [:state, :shock, :data, :param, :total]
            output_var = Symbol(:decomp, comp, class)
            decomp[output_var] = zeros(Ny, length(hs))
        end
    end

    for (i, h) in enumerate(hs)
        # Adjust h for (differences in) conditioning
        h_cond = h - T1_old
        @assert h_cond >= 0

        # Decompose into components
        state_comps, shock_comps = decompose_states_shocks(sys_new, s_new_new_Tmk_Tmk, s_new_new_Tmk_T,
                                                           ϵ_new_new_tgT, classes, T0, k_cond, h_cond)
        check && @assert check_states_shocks_decomp(sys_new, s_new_new_tgt, s_new_new_tgT, classes,
                                                    k_cond, h_cond, state_comps, shock_comps, atol = atol)

        data_comps = decompose_data_revisions(sys_new, s_new_new_Tmk_Tmk, s_old_new_Tmk_Tmk,
                                              classes, T0, k_cond, h_cond)
        param_comps = decompose_param_reest(sys_new, sys_old, s_old_new_Tmk_Tmk, s_old_old_Tmk_Tmk,
                                            classes, T0, k_cond, h_cond)

        for class in classes
            decomp[Symbol(:decompstate, class)][:, i] = state_comps[class]
            decomp[Symbol(:decompshock, class)][:, i] = shock_comps[class]
            decomp[Symbol(:decompdata, class)][:,  i] = data_comps[class]
            decomp[Symbol(:decompparam, class)][:, i] = param_comps[class]

            # Compute total difference
            decomp[Symbol(:decomptotal, class)][:, i] =
                state_comps[class] + shock_comps[class] + data_comps[class] + param_comps[class]
        end
    end

    check && @assert check_total_decomp(m_new, m_old, df_new, df_old, sys_new, sys_old, cond_new, cond_old, classes, decomp;
                                        hs = hs, atol = atol)

    return decomp
end

function decomposition_periods(m_new::M, m_old::M,
                               cond_new::Symbol, cond_old::Symbol) where M<:AbstractModel
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

function prepare_decomposition!(m_new::M, m_old::M, df_new::DataFrame, df_old::DataFrame,
                                params_new::Vector{Float64}, params_old::Vector{Float64},
                                cond_new::Symbol, cond_old::Symbol,
                                k::Int; atol::Float64 = 1e-8) where M<:AbstractModel
    # Update parameters
    DSGE.update!(m_new, params_new)
    DSGE.update!(m_old, params_old)
    sys_new = compute_system(m_new)
    sys_old = compute_system(m_old)

    # Filter and smooth
    s_new_new_tgt = DSGE.filter(m_new, df_new, sys_new, cond_type = cond_new, outputs = [:filt],
                                include_presample = false)[:s_filt]
    s_new_new_tgT, ϵ_new_new_tgT = smooth(m_new, df_new, sys_new, cond_type = cond_new, draw_states = false)

    s_old_new_Tmk_Tmk = DSGE.filter(m_new, df_old, sys_new, cond_type = cond_old, outputs = Symbol[],
                                    include_presample = false)[:s_T]
    s_old_old_Tmk_Tmk = DSGE.filter(m_old, df_old, sys_old, cond_type = cond_old, outputs = Symbol[],
                                    include_presample = false)[:s_T]

    # Check sizes
    T = size(s_new_new_tgt, 2)
    @assert size(s_new_new_tgT, 2) == size(ϵ_new_new_tgT, 2) == T
    @assert isapprox(s_new_new_tgt[:, end], s_new_new_tgT[:, end], atol = atol)

    return sys_new, sys_old, s_new_new_tgt, s_new_new_tgT, ϵ_new_new_tgT, s_old_new_Tmk_Tmk, s_old_old_Tmk_Tmk
end

function decompose_states_shocks(sys_new::System, s_Tmk_Tmk::Vector{Float64},
                                 s_Tmk_T::Vector{Float64}, ϵ_tgT::Matrix{Float64},
                                 classes::Vector{Symbol}, T0::Int, k::Int, h::Int)
    # New parameters, new data
    # State and shock components = y_{T-k+h|T} - y_{T-k+h|T-k}
    #     = Z [ T^h (s_{T-k|T} - s_{T-k|T-k}) + sum_{j=1}^(min(k,h)) ( T^(h-j) R ϵ_{T-k+j|T} ) ]
    #     = state component + shock component
    TTT, RRR, CCC = sys_new[:TTT], sys_new[:RRR], sys_new[:CCC]

    # State component = Z T^h (s_{T-k|T} - s_{T-k|T-k})
    _state_comp = TTT^h * (s_Tmk_T - s_Tmk_Tmk)

    # Shock component = Z sum_{j=1}^(min(k,h)) ( T^(h-j) R ϵ_{T-k+j|T} )
    _shock_comp = zeros(size(s_Tmk_T))
    for j = 1:min(k, h)
        _shock_comp .+= TTT^(h-j) * RRR * ϵ_tgT[:, end-k+j]
    end

    # Map states to desired classes
    state_comps = Dict{Symbol, Vector{Float64}}()
    shock_comps = Dict{Symbol, Vector{Float64}}()
    for class in classes
        ZZ, DD = class_measurement_matrices(sys_new, class)
        state_comps[class] = ZZ * _state_comp
        shock_comps[class] = ZZ * _shock_comp
    end

    return state_comps, shock_comps
end

function decompose_data_revisions(sys_new::System, s_new_Tmk_Tmk::Vector{Float64},
                                  s_old_Tmk_Tmk::Vector{Float64}, classes::Vector{Symbol},
                                  T0::Int, k::Int, h::Int)
    # New parameters, new and old data
    # Data revision component = y^new_{T-k+h|T-k} - y^old_{T-k+h|T-k}
    #     = Z T^h ( s^new_{T-k|T-k} - s^old_{T-k|T-k} )
    TTT = sys_new[:TTT]
    _data_comp = TTT^h * (s_new_Tmk_Tmk - s_old_Tmk_Tmk)

    # Map states to desired classes
    data_comps = Dict{Symbol, Vector{Float64}}()
    for class in classes
        ZZ, _ = class_measurement_matrices(sys_new, class)
        data_comps[class] = ZZ * _data_comp
    end
    return data_comps
end

function decompose_param_reest(sys_new::System, sys_old::System,
                               s_new_Tmk_Tmk::Vector{Float64}, s_old_Tmk_Tmk::Vector{Float64},
                               classes::Vector{Symbol}, T0::Int, k::Int, h::Int)
    # New and old parameters, old data
    # Parameter re-estimation component = y^new_{T-k+h|T-k} - y^old_{T-k+h|T-k}
    TTT_new, CCC_new = sys_new[:TTT], sys_new[:CCC]
    TTT_old, CCC_old = sys_old[:TTT], sys_old[:CCC]

    # y^new_{T-k+h|T-k} = Z^new (T^new^h s^new_{T-k|T-k} + \sum_{j=1}^h (T^new)^(j-1) C^new) + D^new
    s_new_Tmkph_Tmk = TTT_new^h * s_new_Tmk_Tmk
    for j = 1:h
        s_new_Tmkph_Tmk .+= TTT_new^(j-1) * CCC_new
    end

    # y^old_{T-k+h|T-k} = Z^old (T^old^h s^old_{T-k|T-k} + \sum_{j=1}^h (T^old)^(j-1) C^old) + D^old
    s_old_Tmkph_Tmk = TTT_old^h * s_old_Tmk_Tmk
    for j = 1:h
        s_old_Tmkph_Tmk .+= TTT_old^(j-1) * CCC_old
    end

    # Map to desired classes
    param_comps = Dict{Symbol, Vector{Float64}}()
    for class in classes
        ZZ_new, DD_new = class_measurement_matrices(sys_new, class)
        ZZ_old, DD_old = class_measurement_matrices(sys_old, class)

        y_new_Tmkph_Tmk = ZZ_new*s_new_Tmkph_Tmk + DD_new
        y_old_Tmkph_Tmk = ZZ_old*s_old_Tmkph_Tmk + DD_old

        param_comps[class] = y_new_Tmkph_Tmk - y_old_Tmkph_Tmk
    end
    return param_comps
end