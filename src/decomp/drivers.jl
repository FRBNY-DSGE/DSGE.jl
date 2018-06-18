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