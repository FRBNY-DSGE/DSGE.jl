"""
```
decompose_forecast(m_new, m_old, df_new, df_old, input_type, cond_new, cond_old,
    classes; verbose = :low, kwargs...)

decompose_forecast(m_new, m_old, df_new, df_old, params_new, params_old,
    cond_new, cond_old, classes; individual_shocks = false,
    check = false)
```

### Inputs

- `m_new::M` and `m_old::M` where `M<:AbstractModel`
- `df_new::DataFrame` and `df_old::DataFrame`
- `cond_new::Symbol` and `cond_old::Symbol`
- `classes::Vector{Symbol}`: some subset of `[:states, :obs, :pseudo]`

**Method 1 only:**

- `input_type::Symbol`: estimation type to use. Parameters will be loaded using
  `load_draws(m_new, input_type)` and `load_draws(m_old, input_type)` in this
  method

**Method 2 only:**

- `params_new::Vector{Float64}` and `params_old::Vector{Float64}`: single
  parameter draws to use

### Keyword Arguments

- `check::Bool`: whether to check that the individual components add up to the
  correct total difference in forecasts. This roughly doubles the runtime
- `individual_shocks::Bool`: whether to decompose the shock component into
  individual shocks

**Method 1 only:**

- `verbose::Symbol`

### Outputs

The first method returns nothing. The second method returns
`decomp::Dict{Symbol, Matrix{Float64}}`, which has keys of the form
`:decomp<component><class>` and values of size `Ny` x `Nh`, where

- `Ny` is the number of variables in the given class
- `Nh` is the number of common forecast periods, i.e. periods between
  `date_forecast_start(m_new)` and `date_forecast_end(m_old)`
"""
function decompose_forecast(m_new::M, m_old::M, df_new::DataFrame, df_old::DataFrame,
                            input_type::Symbol, cond_new::Symbol, cond_old::Symbol,
                            classes::Vector{Symbol};
                            verbose::Symbol = :low, kwargs...) where M<:AbstractModel
    # Get output file names
    decomp_output_files = get_decomp_output_files(m_new, m_old, input_type, cond_new, cond_old, classes)

    if DSGE.VERBOSITY[verbose] >= DSGE.VERBOSITY[:low]
        info("Decomposing forecast...")
        println("Start time: $(now())")
    end

    # Set up call to lower-level method
    f(params_new::Vector{Float64}, params_old::Vector{Float64}) =
      decompose_forecast(m_new, m_old, df_new, df_old, params_new, params_old,
                         cond_new, cond_old, classes; kwargs...)

    # Single-draw forecasts
    if input_type in [:mode, :mean, :init]

        params_new = load_draws(m_new, input_type, verbose = verbose)
        params_old = load_draws(m_old, input_type, verbose = verbose)
        decomps = f(params_new, params_old)
        write_forecast_decomposition(m_new, m_old, input_type, classes, decomp_output_files, decomps,
                                     verbose = verbose)

    # Multiple-draw forecasts
    elseif input_type == :full

        block_inds, block_inds_thin = forecast_block_inds(m_new, input_type)
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
            decomps = mapfcn(f, params_new, params_old)

            # Assemble outputs from this block and write to file
            decomps = convert(Vector{Dict{Symbol, Array{Float64}}}, decomps)
            decomps = assemble_block_outputs(decomps)
            write_forecast_decomposition(m_new, m_old, input_type, classes, decomp_output_files, decomps,
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
                            individual_shocks::Bool = false, check::Bool = false) where M<:AbstractModel
    # Check numbers of periods
    check_decomposition_periods(m_new, m_old, df_new, df_old, cond_new, cond_old)

    # Forecast
    keep_startdate = date_forecast_start(m_new) # T+1
    keep_enddate   = date_forecast_end(m_old)   # T+H
    H = DSGE.subtract_quarters(keep_enddate, keep_startdate) + 1
    shockdec_splitdate = date_mainsample_end(m_old) # T-k

    f(m::AbstractModel, df::DataFrame, params::Vector{Float64}, cond_type::Symbol; kwargs...) =
        decomposition_forecast(m, df, params, cond_type, keep_startdate, keep_enddate,
                               shockdec_splitdate; kwargs..., check = check)

    out1 = f(m_new, df_new, params_new, cond_new, outputs = [:shockdec])            # new data, new params
    out2 = f(m_old, df_old, params_new, cond_old, outputs = [:shockdec, :forecast]) # old data, new params
    out3 = f(m_old, df_old, params_old, cond_old, outputs = [:forecast])            # old data, old params

    # Initialize output dictionary
    decomp = Dict{Symbol, Array{Float64}}()

    # Decomposition
    for class in classes
        # All elements of out are of size Ny x Nh, where the second dimension
        # ranges from t = T+1:T+H
        forecastvar = Symbol(:forecast, class) # Z s_{T+h} + D
        trendvar    = Symbol(:trend,    class) # Z D
        dettrendvar = Symbol(:dettrend, class) # Z T^{T+h} s_0 + D
        revisionvar = Symbol(:revision, class) # Z \sum_{t=1}^{T-k} T^{T+h-t} R ϵ_t + D
        newsvar     = Symbol(:news,     class) # Z \sum_{t=T-k+1}^{T+h} T^{T+h-t} R ϵ_t + D

        # Data revision
        data_comp = (out1[dettrendvar] - out2[dettrendvar]) + (out1[revisionvar] - out2[revisionvar])
        decomp[Symbol(:decompdata, class)] = data_comp

        # News
        news_comp = out1[newsvar] - out2[newsvar]
        decomp[Symbol(:decompnews, class)] = news_comp

        # Parameter re-estimation
        para_comp = out2[forecastvar] - out3[forecastvar]
        decomp[Symbol(:decomppara, class)] = para_comp

        # Total difference
        total_diff = para_comp + data_comp + news_comp
        decomp[Symbol(:decomptotal, class)] = total_diff
        if check
            @assert total_diff ≈ out1[forecastvar] - out3[forecastvar]
        end
    end

    return decomp
end

"""
```
check_decomposition_periods(m_new, m_old, df_new, df_old, cond_new, cond_old)
```

Check that datasets have the right number of periods.
"""
function check_decomposition_periods(m_new::M, m_old::M, df_new::DataFrame, df_old::DataFrame,
                                     cond_new::Symbol, cond_old::Symbol) where M<:AbstractModel
    # Number of presample periods T0 must be the same
    T0 = n_presample_periods(m_new)
    @assert n_presample_periods(m_old) == T0

    # New model has T main-sample periods
    # Old model has T-k main-sample periods
    T = n_mainsample_periods(m_new)
    k = subtract_quarters(date_forecast_start(m_new), date_forecast_start(m_old))
    @assert k >= 0

    # Number of conditional periods T1 may differ
    T1_new = cond_new == :none ? 0 : n_conditional_periods(m_new)
    T1_old = cond_old == :none ? 0 : n_conditional_periods(m_old)

    # Check DataFrame sizes
    @assert size(df_new, 1) == T0 + T + T1_new
    @assert size(df_old, 1) == T0 + T - k + T1_old

    return nothing
end

"""
```
decomposition_forecast(m, df, params, cond_type, keep_startdate, keep_enddate, shockdec_splitdate;
    outputs = [:forecast, :shockdec], check = false)
```

Equivalent of `forecast_one_draw` for forecast decomposition. `keep_startdate =
date_forecast_start(m_new)` corresponds to time T+1, `keep_enddate =
date_forecast_end(m_old)` to time T+H, and `shockdec_splitdate =
date_mainsample_end(m_old)` to time T-k.

Returns `out::Dict{Symbol, Array{Float64}}`, which has keys determined as follows:

- If `:forecast in outputs` or `check = true`:
  - `:forecast<class>`

- If `:shockdec in outputs`:
  - `:trend<class>`
  - `:dettrend<class>`
  - `:revision<class>`: like a shockdec, but only applying smoothed shocks up to `shockdec_splitdate`
  - `:news<class>`: like a shockdec, but only applying smoothed shocks after `shockdec_splitdate`
"""
function decomposition_forecast(m::AbstractModel, df::DataFrame, params::Vector{Float64}, cond_type::Symbol,
                                keep_startdate::Date, keep_enddate::Date, shockdec_splitdate::Date;
                                outputs::Vector{Symbol} = [:forecast, :shockdec], check::Bool = false)
    # Compute state space
    DSGE.update!(m, params)
    system = compute_system(m)

    # Initialize output dictionary
    out = Dict{Symbol, Array{Float64}}()

    # Smooth
    histstates, histshocks, histpseudo, s_0 = smooth(m, df, system, cond_type = cond_type, draw_states = false)

    # Forecast
    if :forecast in outputs || check
        s_T = histstates[:, end]
        _, forecastobs, forecastpseudo, _ =
            forecast(m, system, s_T, cond_type = cond_type, enforce_zlb = false, draw_shocks = false)

        T = n_mainsample_periods(m)
        if cond_type in [:full, :semi]
            out[:forecastobs]    = DSGE.transplant_forecast_observables(histstates, forecastobs, system, T)
            out[:forecastpseudo] = DSGE.transplant_forecast(histpseudo, forecastpseudo, T)
        else
            out[:forecastobs]    = forecastobs
            out[:forecastpseudo] = forecastpseudo
        end

        # Keep only forecasted periods between keep_startdate and keep_enddate
        forecast_dates = DSGE.quarter_range(date_forecast_start(m), date_forecast_end(m))
        keep_inds = keep_startdate .<= forecast_dates .<= keep_enddate
        out[:forecastobs]    = out[:forecastobs][:, keep_inds]
        out[:forecastpseudo] = out[:forecastpseudo][:, keep_inds]
    end

    # Compute trend, dettrend, and shockdecs
    if :shockdec in outputs
        m <= Setting(:shockdec_startdate, Nullable(keep_startdate))
        m <= Setting(:shockdec_enddate,   Nullable(keep_enddate))

        _, out[:trendobs],    out[:trendpseudo]    = trends(system)
        _, out[:dettrendobs], out[:dettrendpseudo] = deterministic_trends(m, system, s_0)

        # Applying ϵ_{1:T-k}
        t0 = date_mainsample_start(m)
        t1 = shockdec_splitdate
        _, out[:revisionobs], out[:revisionpseudo] =
            DSGE.shock_decompositions(m, system, histshocks, shock_start_date = t0, shock_end_date = t1)

        # Applying ϵ_{T-k+1:end}
        t0 = DSGE.iterate_quarters(shockdec_splitdate, 1)
        t1 = cond_type == :none ? date_mainsample_end(m) : date_conditional_end(m)
        _, out[:newsobs], out[:newspseudo] =
            DSGE.shock_decompositions(m, system, histshocks, shock_start_date = t0, shock_end_date = t1)

        # Shockdec output starts as size Ny x Nh x Ne
        # Sum over shock dimension to get size Ny x Nh
        for output_var in [:revisionobs, :revisionpseudo, :newsobs, :newspseudo]
            out[output_var] = squeeze(sum(out[output_var], 3), 3)
        end
    end

    # Return
    return out
end