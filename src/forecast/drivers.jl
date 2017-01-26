"""
```
forecast_all(m, input_types, cond_types, output_vars; verbose = :low,
    procs = [myid()])
```

Compute forecasts for all specified combinations of conditional data, input
types, and output types.

### Inputs

- `m::AbstractModel`: model object

- `input_types::Vector{Symbol}`: which set of parameters to use, any combination of

```
  - `:mode`: forecast using the modal parameters only
  - `:mean`: forecast using the mean parameters only
  - `:init`: forecast using the initial parameter values only
  - `:full`: forecast using all parameters (full distribution)
  - `:subset`: forecast using a well-defined user-specified subset of draws
```

- `cond_types::Vector{Symbol}`: conditional data types, any combination of

```
  - `:none`: no conditional data
  - `:semi`: use \"semiconditional data\" - average of quarter-to-date
    observations for high frequency series
  - `:full`: use \"conditional data\" - semiconditional plus nowcasts for
    desired observables
```

- `output_vars::Vector{Symbol}`: vector of desired output variables. See
  `forecast_one` for documentation of all possible `output_vars`

### Keyword Arguments

- `verbose::Symbol`: desired frequency of function progress messages printed to
  standard out. One of:

```
  - `:none`: No status updates will be reported.
  - `:low`: Function will report forecast start and input/output locations, as
    well as when each step (preparing inputs, filtering and smoothing,
    forecasting, and shock decompositions) is started.
  - `:high`: Function will report everything in `:low`, as well as elapsed
    times for each step and the file names being written to.
```

- `procs::Vector{Int}`: list of worker processes that have been previously added
  by the user. Defaults to `[myid()]`

### Outputs

None. Output is saved to files returned by
`get_forecast_output_files(m, input_type, cond_type, output_vars)`
for each combination of `input_type` and `cond_type`.
"""
function forecast_all(m::AbstractModel,
                      input_types::Vector{Symbol} = Vector{Symbol}(),
                      cond_types::Vector{Symbol}  = Vector{Symbol}(),
                      output_vars::Vector{Symbol} = Vector{Symbol}();
                      verbose::Symbol             = :low,
                      procs::Vector{Int}          = [myid()])

    for input_type in input_types

        # Load draws and solve system for each draw
        params  = load_draws(m, input_type; verbose = verbose, procs = procs)
        systems = prepare_systems(m, input_type, params; procs = procs)

        # Are we only running IRFs?
        irfs_only = all(x -> contains(string(x), "irf"), output_vars)

        for cond_type in cond_types

            if irfs_only
                # If only running IRFs, don't need to load data or run Kalman filter
                df   = DataFrame()
                kals = dinit(Kalman{S}, 0)
            else
                # Load data
                df = load_data(m; cond_type = cond_type, try_disk = true, verbose = :none)

                # Run Kalman filter to get s_{T|T}
                kals = filter_all(m, df, systems; cond_type = cond_type, procs = procs)
            end

            # Call forecast_one
            forecast_one(m, input_type, cond_type, output_vars;
                         df = df, systems = systems, kals = kals,
                         verbose = verbose, procs = [myid()])
        end
    end
end

"""
`load_draws(m, input_type; subset_inds = 1:0, verbose = :low, procs = [myid()])`

Load and return parameter draws from Metropolis-Hastings.

### Inputs

- `m::AbstractModel`: model object
- `input_type::Symbol`: one of the options for `input_type` described in the
  documentation for `forecast_all`

### Keyword Arguments

- `subset_inds::Range{Int64}`: indices specifying the draws we want to use. See
  `forecast_one` for more detail.
- `verbose::Symbol`: desired frequency of function progress messages printed to
  standard out. One of `:none`, `:low`, or `:high`. If `:low` or greater, prints
  location of input file.
- `procs::Vector{Int}`: list of worker processes over which to distribute
  draws. Defaults to `[myid()]`.

### Outputs

- `params::Matrix{Float64}`: matrix of size `nsim` x `nparams` of parameter
  draws

where `nsim` is the number of draws saved in Metropolis-Hastings.

### Notes

If `nsim` is not divisible by `jstep * nprocs`, where:

- `jstep = get_jstep(m, nsim)` is the thinning step size for the forecast step
- `nprocs = length(procs)` is the number of processes over which we distribute draws

then we truncate the draws so that `mod(nsim_new, jstep * nprocs) == 0`.
"""
function load_draws(m::AbstractModel, input_type::Symbol;
    subset_inds::Range{Int64} = 1:0,
    verbose::Symbol = :low, procs::Vector{Int} = [myid()])

    input_file_name = get_input_file(m, input_type)
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("Loading draws from $input_file_name")
    end

    # Read infiles
    if input_type in [:mean, :mode]
        tmp = map(Float64, h5read(input_file_name, "params"))
        params = reshape(tmp, 1, size(tmp, 1))

    elseif input_type in [:full, :subset, :block]
        if input_type == :full
            params = map(Float64, h5read(input_file_name, "mhparams"))
        else
            if isempty(subset_inds)
                error("Must supply nonempty range of subset_inds if input_type = subset or block")
            else
                params = map(Float64, h5read(input_file_name, "mhparams", (subset_inds, :)))
            end
        end

        # Truncate number of draws if necessary
        nsim = size(params,1)
        jstep = get_jstep(m, nsim)
        nprocs = length(procs)
        remainder = mod(nsim, jstep * nprocs)
        if remainder != 0
            nsim_new = nsim - remainder
            warn("Number of draws read in, $nsim, is not divisible by jstep * nprocs = $(jstep * nprocs). Taking the first $nsim_new draws instead.")
            params = params[1:nsim_new, :]
        end

    elseif input_type in [:init]
        init_parameters!(m)
        tmp = Float64[α.value for α in m.parameters]
        params = reshape(tmp, 1, size(tmp,1))
    end

    return params
end

"""
```
prepare_systems(m, input_type, params; procs = [myid()])
```

Returns a `DVector{System{Float64}}` of `System` objects, one for each draw in
`params` (after thinning), by recomputing the state-space system.

### Inputs

- `m::AbstractModel`: model object
- `input_type`: one of the options for `input_type` described in the
  documentation for `forecast_all`.
- `params::Matrix{Float64}`: matrix of size `nsim` x `nparams` of parameter
  draws

### Keyword Arguments

- `procs::Vector{Int}`: list of worker processes over which to distribute
  draws. Defaults to `[myid()]`.

### Outputs

- `systems::DVector{System{Float64}}`: vector of `n_sim_forecast` many `System`
  objects, one for each draw

where `n_sim_forecast = n_sim / jstep` is the number of draws after thinning a
second time.
"""
function prepare_systems(m::AbstractModel, input_type::Symbol,
    params::Matrix{Float64}; procs::Vector{Int} = [myid()])

    # Reset procs to [myid()] if necessary
    procs = reset_procs(m, procs, Nullable(input_type))

    # Setup
    n_sim = size(params,1)
    jstep = get_jstep(m, n_sim)
    n_sim_forecast = convert(Int, n_sim/jstep)

    if input_type in [:mean, :mode, :init]
        update!(m, vec(params))
        systems = dfill(compute_system(m), (1,), procs)
    elseif input_type in [:full, :subset, :block]
        nprocs = length(procs)
        systems = DArray((n_sim_forecast,), procs, [nprocs]) do I
            draw_inds = first(I)
            ndraws_local = convert(Int, n_sim_forecast/nprocs)
            localpart = Vector{System{Float64}}(ndraws_local)

            for i in draw_inds
                j = i * jstep
                i_local = mod(i-1, ndraws_local) + 1

                params_j = vec(params[j,:])
                update!(m, params_j)
                localpart[i_local] = compute_system(m)
            end
            return localpart
        end
    else
        throw(ArgumentError("Not implemented."))
    end

    return systems
end

"""
```
forecast_one(m, input_type, cond_type, output_vars;
    df = DataFrame(), systems = dinit(System{Float64}, 0},
    kals = dinit(Kalman{Float64}, 0), subset_inds = 1:0,
    forecast_string = "", verbose = :low, procs = [myid()])
```

Compute, save, and return `output_vars` for input draws given by `input_type`
and conditional data case given by `cond_type`.

### Inputs

- `m::AbstractModel`: model object
- `input_type::Symbol`: See documentation for `forecast_all`
- `cond_type::Symbol`: See documentation for `forecast_all`
- `output_vars::Vector{Symbol}`: vector of desired output variables. See Outputs
  section

### Keyword Arguments

- `df::DataFrame`: Historical data. If `cond_type in [:semi, :full]`, then the
   final row of `df` should be the period containing conditional data. If not
   provided, will be loaded using `load_data` with the appropriate `cond_type`
- `systems::DVector{System{Float64}}`: vector of `n_sim_forecast` many `System`
  objects, one for each draw. If not provided, will be loaded using
  `prepare_systems`
- `kals::DVector{Kalman{Float64}}`: vector of `n_sim_forecast` many `Kalman`
  objects. If not provided, will be loaded using `filter_all`
- `subset_inds::Range{Int64}`: indices specifying the draws we want to use. If a
  more sophisticated selection criterion is desired, the user is responsible for
  determining the indices corresponding to that criterion. If `input_type` is
  not `subset`, `subset_inds` will be ignored
- `forecast_string::AbstractString`: short string identifying the subset to be
  appended to the output filenames. If `input_type = :subset` and
  `forecast_string` is empty, an error is thrown.
- `verbose::Symbol`: desired frequency of function progress messages printed to
  standard out. One of:

```
  - `:none`: No status updates will be reported.
  - `:low`: Function will report forecast start and input/output locations, as
    well as when each step (preparing inputs, filtering and smoothing,
    forecasting, and shock decompositions) is started.
  - `:high`: Function will report everything in `:low`, as well as elapsed
    times for each step and the file names being written to.
```

- `procs::Vector{Int}`: list of worker processes that have been previously added
  by the user. Defaults to `[myid()]`

### Output

- `forecast_outputs::Dict{Symbol, DArray{Float64}}`: dictionary of forecast
  outputs. Keys are `output_vars`, which is some subset of:

```
  - `:histstates`: `DArray{Float64, 3}` of smoothed historical states
  - `:histobs`: `DArray{Float64, 3}` of smoothed historical data
  - `:histpseudo`: `DArray{Float64, 3}` of smoothed historical
    pseudo-observables (if a pseudo-measurement equation has been provided for
    this model type)
  - `:histshocks`: `DArray{Float64, 3}` of smoothed historical shocks
  - `:forecaststates`: `DArray{Float64, 3}` of forecasted states
  - `:forecastobs`: `DArray{Float64, 3}` of forecasted observables
  - `:forecastpseudo`: `DArray{Float64, 3}` of forecasted pseudo-observables (if
    a pseudo-measurement equation has been provided for this model type)
  - `:forecastshocks`: `DArray{Float64, 3}` of forecasted shocks
  - `:bddforecaststates`, `:bddforecastobs`, `:bddforecastpseudo`, and
    `:bddforecastshocks`: `DArray{Float64, 3}`s of forecasts where we enforce
    the zero lower bound to be `forecast_zlb_value(m)`
  - `:shockdecstates`: `DArray{Float64, 4}` of state shock decompositions
  - `:shockdecobs`: `DArray{Float64, 4}` of observable shock decompositions
  - `:shockdecpseudo`: `DArray{Float64, 4}` of pseudo-observable shock
    decompositions (if a pseudo-measurement equation has been provided for this
    model type)
  - `:irfstates`: `DArray{Float64, 4}` of state impulse responses
  - `:irfobs`: `DArray{Float64, 4}` of observable impulse responses
  - `:irfpseudo`: `DArray{Float64, 4}` of pseudo observable impulse responses
      (if a pseudo-measurement equation has been provided)
```

### Notes

- `forecast_one` prepares inputs to the forecast for a particular combination of
  input, output and cond types by distributing system matrices and final
  historical state vectors across `procs`. The user must add processes
  separately and pass the worker identification numbers here.
- if `ndraws % length(procs) != 0`, `forecast_one` truncates `procs` so that the
  draws can be evenly distributed across workers. This is required by the
  `DistributedArrays` package.
"""
function forecast_one(m::AbstractModel{Float64},
    input_type::Symbol, cond_type::Symbol, output_vars::Vector{Symbol};
    df::DataFrame = DataFrame(),
    systems::DVector{System{Float64}} = dinit(System{Float64}, 0),
    kals::DVector{Kalman{Float64}} = dinit(Kalman{Float64}, 0),
    block_number::Nullable{Int64} = Nullable{Int64}(),
    subset_inds::Range{Int64} = 1:0,
    forecast_string::AbstractString = "", verbose::Symbol = :low,
    procs::Vector{Int} = [myid()])

    ### 0. Setup

    # Compute everything that will be needed to plot original output_vars
    output_vars = add_requisite_output_vars(output_vars)

    # Prepare forecast outputs
    forecast_output = Dict{Symbol, DArray{Float64}}()
    forecast_output_files = get_forecast_output_files(m, input_type, cond_type, output_vars;
                                forecast_string = forecast_string)
    output_dir = rawpath(m, "forecast")

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println()
        if input_type == :block
            block_num = get(block_number)
            info("Forecasting block $(block_num) of $(n_forecast_blocks(m))...")
        else
            info("Forecasting input_type = $input_type, cond_type = $cond_type...")
            println("Start time: $(now())")
            println("Forecast outputs will be saved in $output_dir")
        end
    end

    # If forecasting in blocks, call forecast_one with input_type = :block
    if forecast_blocking(m) && input_type == :full
        block_inds = forecast_block_inds(m; procs = procs)
        nblocks = n_forecast_blocks(m)
        block_verbose = verbose == :none ? :none : :low
        total_forecast_time = 0.0

        for block = 1:nblocks
            tic()
            forecast_one(m, :block, cond_type, output_vars;
                df = df, block_number = Nullable(block),
                subset_inds = block_inds[block],
                verbose = block_verbose, procs = procs)
            darray_closeall()
            gc()

            # Calculate time to complete this block, average block time, and
            # expected time to completion
            if VERBOSITY[verbose] >= VERBOSITY[:low]
                block_time = toq()
                total_forecast_time += block_time
                total_forecast_time_min = total_forecast_time/60
                expected_time_remaining     = (total_forecast_time/block)*(nblocks - block)
                expected_time_remaining_min = expected_time_remaining/60

                println("\nCompleted $block of $nblocks blocks.")
                println("Total time to compute $block blocks: $total_forecast_time_min minutes")
                println("Expected time remaining in forecast: $expected_time_remaining_min minutes")
            end
        end
        return forecast_output
    end

    # Check that provided forecast inputs are well-formed
    # If an input is not provided, load it
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("\nPreparing forecast inputs...")
    end
    @time_verbose df, systems, kals, procs =
        prepare_forecast_inputs!(m, input_type, cond_type, output_vars; df = df,
            systems = systems, kals = kals, subset_inds = subset_inds,
            verbose = verbose, procs = procs)

    nprocs = length(procs)
    ndraws = length(systems)


    ### 1. Smoothed Histories

    # Must run smoother for conditional data in addition to explicit cases
    hist_vars = [:histstates, :histpseudo, :histshocks]
    shockdec_vars = [:shockdecstates, :shockdecpseudo, :shockdecobs]
    dettrend_vars = [:dettrendstates, :dettrendpseudo, :dettrendobs]
    smooth_vars = vcat(hist_vars, shockdec_vars, dettrend_vars)

    hists_to_compute = intersect(output_vars, hist_vars)
    if !isempty(intersect(output_vars, smooth_vars)) || cond_type in [:semi, :full]

        if VERBOSITY[verbose] >= VERBOSITY[:low]
            println("\nSmoothing $(hists_to_compute)...")
        end
        @time_verbose histstates, histshocks, histpseudo, initial_states =
            smooth_all(m, df, systems, kals; cond_type = cond_type, procs = procs)

        # For conditional data, transplant the obs/state/pseudo vectors from hist to forecast
        if cond_type in [:full, :semi]
            T = n_mainsample_periods(m)

            forecast_output[:histstates] = transplant_history(histstates, T)
            forecast_output[:histshocks] = transplant_history(histshocks, T)
	        if :histpseudo in output_vars
                forecast_output[:histpseudo] = transplant_history(histpseudo, T)
            end
        else
            forecast_output[:histstates] = histstates
            forecast_output[:histshocks] = histshocks
	        if :histpseudo in output_vars
                forecast_output[:histpseudo] = histpseudo
            end
        end

        write_forecast_outputs(m, hists_to_compute, forecast_output_files, forecast_output;
                               block_number = block_number, verbose = verbose)
    end


    ### 2. Forecasts

    # 2A. Unbounded forecasts

    forecast_vars = [:forecaststates, :forecastobs, :forecastpseudo, :forecastshocks]
    forecasts_to_compute = intersect(output_vars, forecast_vars)
    if !isempty(forecasts_to_compute)
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            println("\nForecasting $(forecasts_to_compute)...")
        end
        @time_verbose forecaststates, forecastobs, forecastpseudo, forecastshocks =
            forecast(m, systems, kals; cond_type = cond_type, enforce_zlb = false,
                     procs = procs)

        # For conditional data, transplant the obs/state/pseudo vectors from hist to forecast
        if cond_type in [:full, :semi]
            forecast_output[:forecaststates] = transplant_forecast(histstates, forecaststates, T)
            forecast_output[:forecastshocks] = transplant_forecast(histshocks, forecastshocks, T)
	        if :forecastpseudo in output_vars
                forecast_output[:forecastpseudo] = transplant_forecast(histpseudo, forecastpseudo, T)
            end
            forecast_output[:forecastobs] = transplant_forecast_observables(histstates, forecastobs, systems, T)
        else
            forecast_output[:forecaststates] = forecaststates
            forecast_output[:forecastshocks] = forecastshocks
	        if :forecastpseudo in output_vars
                forecast_output[:forecastpseudo] = forecastpseudo
            end
            forecast_output[:forecastobs] = forecastobs
        end

        write_forecast_outputs(m, forecasts_to_compute, forecast_output_files, forecast_output;
                               block_number = block_number, verbose = verbose)
    end


    # 2B. Bounded forecasts

    forecast_vars_bdd = [:bddforecaststates, :bddforecastobs, :bddforecastpseudo, :bddforecastshocks]
    forecasts_to_compute = intersect(output_vars, forecast_vars_bdd)
    if !isempty(forecasts_to_compute)
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            println("\nForecasting $(forecasts_to_compute)...")
        end
        @time_verbose forecaststates, forecastobs, forecastpseudo, forecastshocks =
            forecast(m, systems, kals; cond_type = cond_type, enforce_zlb = true,
                     procs = procs)

        # For conditional data, transplant the obs/state/pseudo vectors from hist to forecast
        if cond_type in [:full, :semi]
            forecast_output[:bddforecaststates] = transplant_forecast(histstates, forecaststates, T)
            forecast_output[:bddforecastshocks] = transplant_forecast(histshocks, forecastshocks, T)
	        if :bddforecastpseudo in output_vars
                forecast_output[:bddforecastpseudo] = transplant_forecast(histpseudo, forecastpseudo, T)
            end
            forecast_output[:bddforecastobs] = transplant_forecast_observables(histstates, forecastobs, systems, T)
        else
            forecast_output[:bddforecaststates] = forecaststates
            forecast_output[:bddforecastshocks] = forecastshocks
	        if :bddforecastpseudo in output_vars
                forecast_output[:bddforecastpseudo] = forecastpseudo
            end
            forecast_output[:bddforecastobs] = forecastobs
        end

        write_forecast_outputs(m, forecasts_to_compute, forecast_output_files, forecast_output;
                               block_number = block_number, verbose = verbose)
    end


    ### 3. Shock Decompositions

    shockdecs_to_compute = intersect(output_vars, shockdec_vars)
    if !isempty(shockdecs_to_compute)
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            println("\nComputing shock decompositions for $(shockdecs_to_compute)...")
        end
        @time_verbose shockdecstates, shockdecobs, shockdecpseudo =
            shock_decompositions(m, systems, histshocks; procs = procs)

        forecast_output[:shockdecstates] = shockdecstates
        forecast_output[:shockdecobs]    = shockdecobs
        if :shockdecpseudo in output_vars
            forecast_output[:shockdecpseudo] = shockdecpseudo
        end

        write_forecast_outputs(m, shockdecs_to_compute, forecast_output_files, forecast_output;
                               block_number = block_number, verbose = verbose)
    end


    ### 4. Trend

    trend_vars = vcat([:trendstates, :trendobs, :trendpseudo])
    trends_to_compute = intersect(output_vars, trend_vars)

    if !isempty(trends_to_compute)
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            println("\nComputing trend for $(trends_to_compute)...")
        end
        @time_verbose trendstates, trendobs, trendpseudo =
            trends(m, systems; procs = procs)

        forecast_output[:trendstates] = trendstates
        forecast_output[:trendobs]    = trendobs
        if :trendpseudo in output_vars
            forecast_output[:trendpseudo] = trendpseudo
        end

        write_forecast_outputs(m, trends_to_compute, forecast_output_files, forecast_output;
                               block_number = block_number, verbose = verbose)
    end

    ### 5. Deterministic Trend

    dettrends_to_compute = intersect(output_vars, dettrend_vars)
    if !isempty(dettrends_to_compute)

        # Compute deterministic trend
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            println("\nComputing deterministic trend for $(dettrends_to_compute)...")
        end

        @time_verbose dettrendstates, dettrendobs, dettrendpseudo =
            deterministic_trends(m, systems, initial_states; procs = procs)

        forecast_output[:dettrendstates] = dettrendstates
        forecast_output[:dettrendobs]    = dettrendobs
        if :dettrendpseudo in output_vars
            forecast_output[:dettrendpseudo] = dettrendpseudo
        end

        dettrend_write_vars = [symbol("dettrend$c") for c in unique(map(get_class, dettrends_to_compute))]
        write_forecast_outputs(m, dettrends_to_compute, forecast_output_files, forecast_output;
                               block_number = block_number, verbose = verbose)
    end

    ### 6. Impulse Responses
    irf_vars = [:irfstates, :irfobs, :irfpseudo]
    irfs_to_compute = intersect(output_vars, irf_vars)
    if !isempty(irfs_to_compute)
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            println("\nComputing impulse responses for $(irfs_to_compute)...")
        end
        @time_verbose irfstates, irfobs, irfpseudo = impulse_responses(m, systems; procs = procs)

        forecast_output[:irfstates] = irfstates
        forecast_output[:irfobs] = irfobs
        if :irfpseudo in output_vars
            forecast_output[:irfpseudo] = irfpseudo
        end

        write_forecast_outputs(m, irfs_to_compute, forecast_output_files, forecast_output;
                               block_number = block_number, verbose = verbose)
    end


    # Return only desired output_vars
    for key in keys(forecast_output)
        if !(key in output_vars)
            close(forecast_output[key])   # Explicitly garbage-collect DArray
            delete!(forecast_output, key) # Remove from Dict
        end
    end

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("\nForecast complete: $(now())")
    end

    return forecast_output
end
