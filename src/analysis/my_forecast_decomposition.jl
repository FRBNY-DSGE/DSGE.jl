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

        decomps = decompose_forecast(m_new, m_old, df_new, df_old, params_new, params_old, classes;
                                     hs = hs, cond_new = cond_new, cond_old = cond_old, kwargs...)
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
                             decompose_forecast(m_new, m_old, df_new, df_old, param_new, param_old, classes;
                                                hs = hs, cond_new = cond_new, cond_old = cond_old, kwargs...),
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
                            classes::Vector{Symbol}; hs = 1:forecast_horizons(m_old),
                            cond_new::Symbol = :none, cond_old::Symbol = cond_new,
                            check::Bool = false, atol::Float64 = 1e-8)

    # Compute numbers of periods
    # h_cond and k_cond are h and k adjusted for (differences in) conditioning
    T0, T, k, T1_new, T1_old, k_cond =
        decomposition_periods(m_new, m_old, cond_new = cond_new, cond_old = cond_old)

    if check
        @assert size(df_new, 1) == T0 + T + T1_new
        @assert size(df_old, 1) == T0 + T + T1_old - k
    end

    # Update parameters, filter, and smooth
    sys_new, sys_old, kal_new_new, kal_new_old, kal_old_old, s_tgT, ϵ_tgT =
        prepare_decomposition!(m_new, m_old, df_new, df_old, params_new, params_old,
                               cond_new = cond_new, cond_old = cond_old)

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
            # y^{new,new}_{t|T-k}, t = 1:T:H
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
                # Total difference = y^{new,new}_{T+h-k|T} - y^{old,old}_{T-k+h|T-k}
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

function decomposition_periods(m_new::AbstractModel, m_old::AbstractModel;
                               cond_new::Symbol = :none, cond_old::Symbol = cond_old)
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
                                params_new::Vector{Float64}, params_old::Vector{Float64};
                                cond_new::Symbol = :none, cond_old::Symbol = cond_old)
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
    s_tgT, ϵ_tgT = smooth(m_new, df_new, sys_new, kal_new_new,
                          cond_type = cond_new, draw_states = false)

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
    s_tgt = kal_new[:filt][:, (T0+1):end] # s_{t|t}

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

    s_new_tgt = kal_new[:filt][:, (T0+1):end] # s^{new}_{t|t}
    s_old_tgt = kal_old[:filt][:, (T0+1):end] # s^{old}_{t|t}
    try
        check && @assert size(s_new_tgt, 2) == size(s_old_tgt, 2) + k
    catch ex
        @show size(s_new_tgt, 2)
        @show size(s_old_tgt, 2)
        @show k
        throw(ex)
    end

    # Return
    data_comp = ZZ * TTT^h * (s_new_tgt[:, end-k] - s_old_tgt[:, end-k])
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

    s_new_tgt = kal_new[:filt][:, (T0+1):end] # s^{new}_{t|t}
    s_old_tgt = kal_old[:filt][:, (T0+1):end] # s^{old}_{t|t}
    check && @assert size(s_new_tgt, 2) == size(s_old_tgt, 2)

    # y^{new}_{T-k+h|T-k} = Z^{new} (T^{new}^h s^{new}_{T-k|T-k} + C^{new}) + D^{new}
    s_new_Tmk_Tmk = s_new_tgt[:, end-k]
    y_new_Tmkph_Tmk = ZZ_new * (TTT_new^h * s_new_Tmk_Tmk + CCC_new) + DD_new

    # y^{old}_{T-k+h|T-k} = Z^{old} (T^{old}^h s^{old}_{T-k|T-k} + C^{old}) + D^{old}
    s_old_Tmk_Tmk = s_old_tgt[:, end]
    y_old_Tmkph_Tmk = ZZ_old * (TTT_old^h * s_old_Tmk_Tmk + CCC_old) + DD_old

    # Return
    param_comp = y_new_Tmkph_Tmk - y_old_Tmkph_Tmk
    return param_comp
end



### HELPER FUNCTIONS

function class_system_matrices(sys::System, class::Symbol)
    TTT, RRR, CCC = sys[:TTT], sys[:RRR], sys[:CCC]
    ZZ, DD = if class == :obs
        sys[:ZZ], sys[:DD]
    elseif class == :pseudo
        sys[:ZZ_pseudo], sys[:DD_pseudo]
    elseif class == :states
        I, zeros(size(TTT, 1))
    else
        error("Invalid class: $class. Must be :obs, :pseudo, or :state")
    end
    return ZZ, DD, TTT, RRR, CCC
end

function compute_history_and_forecast(m::AbstractModel, df::DataFrame, class::Symbol;
                                      cond_type::Symbol = :none)
    sys = compute_system(m)
    ZZ, DD, _, _, _ = class_system_matrices(sys, class)

    kal = DSGE.filter(m, df, sys, cond_type = cond_type)
    histstates, _ = smooth(m, df, sys, kal, cond_type = cond_type, draw_states = false)
    hist = ZZ * histstates .+ DD

    fcast = Dict{Symbol, Matrix{Float64}}()
    fcast[:states], fcast[:obs], fcast[:pseudo] = forecast(m, sys, kal[:zend])

    return hcat(hist, fcast[class])
end

function get_decomp_output_files(m_new::AbstractModel, m_old::AbstractModel,
                                 input_type::Symbol, cond_new::Symbol, cond_old::Symbol,
                                 classes::Vector{Symbol})
    output_files = Dict{Symbol, String}()
    for comp in [:state, :shock, :data, :param, :total]
        for class in classes
            output_var = Symbol(:decomp, comp, class)

            fn_new = get_forecast_filename(m_new, input_type, cond_new, output_var)
            fn_old = get_forecast_filename(m_old, input_type, cond_old, output_var)

            dir = dirname(fn_new)
            base_new, ext = splitext(basename(fn_new))
            base_old, _ = splitext(basename(fn_old))
            base_old = replace(base_old, string(output_var), "")
            base_old = spec(m_old) * "_" * subspec(m_old) * base_old

            output_files[output_var] = joinpath(dir, base_new * "__" * base_old * ext)
        end
    end
    return output_files
end

function write_forecast_decomposition(m_new::AbstractModel, m_old::AbstractModel,
                                      input_type::Symbol, classes::Vector{Symbol},
                                      hs::Union{Int, UnitRange{Int}},
                                      decomp_output_files::Dict{Symbol, String},
                                      decomps::Dict{Symbol, Array{Float64}};
                                      block_number::Nullable{Int} = Nullable{Int}(),
                                      block_inds::Range{Int} = 1:0,
                                      verbose::Symbol = :low)
    for comp in [:state, :shock, :data, :param, :total]
        for class in classes
            prod = Symbol(:decomp, comp)
            var = Symbol(prod, class)
            filepath = decomp_output_files[var]

            if isnull(block_number) || get(block_number) == 1
                jldopen(filepath, "w") do file
                    # Write metadata
                    # Pass in m_old because its date_mainsample_end is used to calculate dates
                    DSGE.write_forecast_metadata(m_old, file, prod, class, hs = hs)

                    # Pre-allocate HDF5 dataset which will contain all draws
                    if !isnull(block_number) && get(block_number) == 1
                        # Determine forecast output size
                        dims = DSGE.get_forecast_output_dims(m_new, input_type, Symbol(:forecast, class))
                        dims = (dims[1], dims[2], length(hs)) # dims is ndraws x nvars x nperiods
                        block_size = forecast_block_size(m_new)
                        chunk_dims = collect(dims)
                        chunk_dims[1] = block_size

                        # Initialize dataset
                        pfile = file.plain
                        HDF5.d_create(pfile, "arr", datatype(Float64), dataspace(dims...), "chunk", chunk_dims)
                    end
                end
            end

            jldopen(filepath, "r+") do file
                if isnull(block_number)
                    write(file, "arr", decomps[var])
                else
                    DSGE.write_forecast_block(file, decomps[var], block_inds)
                end
            end

            if DSGE.VERBOSITY[verbose] >= DSGE.VERBOSITY[:high]
                println(" * Wrote $(basename(filepath))")
            end
        end # of loop over classes
    end # of loop over comps
end