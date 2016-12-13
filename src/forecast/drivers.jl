"""
```
forecast_all(m, cond_types, input_types, output_vars; verbose = :low,
    procs = [myid()])
```

Compute forecasts for all specified combinations of conditional data, input
types, and output types.

### Inputs

- `m::AbstractModel`: model object

- `cond_types::Vector{Symbol}`: conditional data types, any combination of

```
  - `:none`: no conditional data
  - `:semi`: use \"semiconditional data\" - average of quarter-to-date
    observations for high frequency series
  - `:full`: use \"conditional data\" - semiconditional plus nowcasts for
    desired observables
```

- `input_types::Vector{Symbol}`: which set of parameters to use, any combination of

```
  - `:mode`: forecast using the modal parameters only
  - `:mean`: forecast using the mean parameters only
  - `:init`: forecast using the initial parameter values only
  - `:full`: forecast using all parameters (full distribution)
  - `:subset`: forecast using a well-defined user-specified subset of draws
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
`get_output_files(m, input_type, output_vars, cond_type)`
for each combination of `input_type` and `cond_type`.
"""
function forecast_all(m::AbstractModel,
                      cond_types::Vector{Symbol}  = Vector{Symbol}(),
                      input_types::Vector{Symbol} = Vector{Symbol}(),
                      output_vars::Vector{Symbol} = Vector{Symbol}();
                      verbose::Symbol             = :low,
                      procs::Vector{Int}          = [myid()])

    for cond_type in cond_types
        df = load_data(m; cond_type=cond_type, try_disk=true, verbose=:none)
        for input_type in input_types
            my_procs = if input_type in [:init, :mode, :mean]
                [myid()]
            else
                procs
            end
            forecast_one(m, df; cond_type=cond_type, input_type=input_type,
                output_vars=output_vars, verbose=verbose, procs=my_procs)
        end
    end
end

"""
`load_draws(m, input_type; subset_inds = [], verbose = :low, procs = [myid()])`

Load and return parameter draws, transition matrices, and final state vectors
from Metropolis-Hastings. Single draws are reshaped to have additional singleton
dimensions, and missing variables without sufficient information are initialized
to null values of appropriate types.

### Inputs

- `m::AbstractModel`: model object
- `input_type::Symbol`: one of the options for `input_type` described in the
  documentation for `forecast_all`

### Keyword Arguments

- `subset_inds::Vector{Int}`: indices specifying the draws we want to use. See
  `forecast_one` for more detail.
- `verbose::Symbol`: desired frequency of function progress messages printed to
  standard out. One of `:none`, `:low`, or `:high`. If `:low` or greater, prints
  location of input file.
- `procs::Vector{Int}`: list of worker processes over which to distribute
  draws. Defaults to `[myid()]`.

### Outputs

- `params::Matrix{Float64}`: matrix of size `nsim` x `nparams` of parameter
  draws
- `TTT::Array{Float64, 3}`: array of size `nsim` x `nstates` x `nstates` of
  transition matrix draws
- `RRR::Array{Float64, 3}`: array of size `nsim` x `nstates` x `nshocks` of
  transition matrix draws
- `CCC::Array{Float64, 3}`: array of size `nsim` x `nstates` x 1 of transition
  matrix draws

where `nsim` is the number of draws saved in Metropolis-Hastings.

### Notes

If `nsim` is not divisible by `jstep * nprocs`, where:

- `jstep = get_jstep(m, nsim)` is the thinning step size for the forecast step
- `nprocs = length(procs)` is the number of processes over which we distribute draws

then we truncate the draws so that `mod(nsim_new, jstep * nprocs) == 0`.
"""
function load_draws(m::AbstractModel, input_type::Symbol;
    subset_inds::Vector{Int} = Vector{Int}(), verbose::Symbol = :low,
    procs::Vector{Int} = [myid()])

    input_file_name = get_input_file(m, input_type)
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("Loading draws from $input_file_name")
    end

    # Read infiles
    if input_type in [:mean, :mode]
        tmp = h5open(input_file_name, "r") do f
            map(Float64, read(f, "params"))
        end
        params = reshape(tmp, 1, size(tmp,1))
        TTT  = Array{Float64}(0,0,0)
        RRR  = Array{Float64}(0,0,0)
        CCC  = Array{Float64}(0,0,0)

    elseif input_type in [:full, :subset]
        if input_type == :full
            h5open(input_file_name, "r") do f
                params = map(Float64, read(f, "mhparams"))
                TTT    = map(Float64, read(f, "mhTTT"))
                RRR    = map(Float64, read(f, "mhRRR"))
                if "mhCCC" in names(f)
                    CCC = map(Float64, read(f, "mhCCC"))
                else
                    CCC = Array{Float64}(0,0,0)
                end
            end
        else
            if isempty(subset_inds)
                error("Must supply nonempty vector of subset_inds if input_type = subset")
            end
            h5open(input_file_name, "r") do f
                params = map(Float64, read(f, "mhparams")[subset_inds, :])
                TTT    = map(Float64, read(f, "mhTTT")[subset_inds, :, :])
                RRR    = map(Float64, read(f, "mhRRR")[subset_inds, :, :])
                if "mhCCC" in names(f)
                    CCC = map(Float64, read(f, "mhCCC")[subset_inds, :])
                else
                    CCC = Array{Float64}(0,0,0)
                end
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
            TTT    = TTT[1:nsim_new, :, :]
            RRR    = RRR[1:nsim_new, :, :]
            if !isempty(CCC)
                CCC = CCC[1:nsim_new, :, :]
            end
        end

    elseif input_type in [:init]
        init_parameters!(m)
        tmp = Float64[α.value for α in m.parameters]
        params = reshape(tmp, 1, size(tmp,1))
        TTT  = Array{Float64}(0,0,0)
        RRR  = Array{Float64}(0,0,0)
        CCC  = Array{Float64}(0,0,0)
    end

    return params, TTT, RRR, CCC
end

"""
```
prepare_systems(m, input_type, params, TTT, RRR, CCC; procs = [myid()])
```

Returns a `DVector{System{Float64}}` of `System` objects constructed from the
given sampling outputs that is suitable for input to forecasting. In the
one-draw case (`input_type in [:mode, :mean, :init]`), the model is re-solved
and the state-space system is recomputed. In the many-draw case (`input_type in
[:full, :subset]`), the outputs from sampling are simply repackaged to the
appropriate shape.

### Inputs

- `m::AbstractModel`: model object
- `input_type`: one of the options for `input_type` described in the
  documentation for `forecast_all`.
- `params::Matrix{Float64}`: matrix of size `nsim` x `nparams` of parameter
  draws
- `TTT::Array{Float64, 3}`: array of size `nsim` x `nstates` x `nstates` of
  transition matrix draws
- `RRR::Array{Float64, 3}`: array of size `nsim` x `nstates` x `nshocks` of
  transition matrix draws
- `CCC::Array{Float64, 3}`: array of size `nsim` x `nstates` x 1 of transition
  matrix draws

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
    params::Matrix{Float64}, TTT::Array{Float64, 3}, RRR::Array{Float64, 3},
    CCC::Array{Float64, 3}; procs::Vector{Int} = [myid()])

    # Reset procs to [myid()] if necessary
    procs = reset_procs(m, procs, Nullable(input_type))

    # Setup
    n_sim = size(params,1)
    jstep = get_jstep(m, n_sim)
    n_sim_forecast = convert(Int, n_sim/jstep)

    if input_type in [:mean, :mode, :init]
        update!(m, vec(params))
        systems = dfill(compute_system(m), (1,), procs)
    elseif input_type in [:full, :subset]
        empty = isempty(CCC)
        nprocs = length(procs);
        systems = DArray((n_sim_forecast,), procs, [nprocs]) do I
            draw_inds = first(I)
            ndraws_local = convert(Int, n_sim_forecast/nprocs)
            localpart = Vector{System{Float64}}(ndraws_local)

            for i in draw_inds
                j = i * jstep
                i_local = mod(i-1, ndraws_local) + 1

                # Prepare transition eq
                TTT_j = squeeze(TTT[j, :, :], 1)
                RRR_j = squeeze(RRR[j, :, :], 1)

                if empty
                    trans_j = Transition(TTT_j, RRR_j)
                else
                    CCC_j = squeeze(CCC[j, :, :], 1)
                    trans_j = Transition(TTT_j, RRR_j, CCC_j)
                end

                # Prepare measurement eq
                params_j = vec(params[j,:])
                update!(m, params_j)
                meas_j   = measurement(m, trans_j; shocks = true)

                # Prepare pseudo-measurement eq
                pseudo_meas_j = if forecast_pseudoobservables(m)
                    _, pseudo_mapping = pseudo_measurement(m)
                    Nullable(pseudo_mapping)
                else
                    Nullable{PseudoObservableMapping{Float64}}()
                end

                # Prepare system
                localpart[i_local] = System(trans_j, meas_j, pseudo_meas_j)
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
prepare_forecast_inputs(m, df; input_type = :mode, cond_type = :none,
    subset_inds = [], verbose = :low, procs = [myid()])
```

Load draws for this input type, prepare a System object for each draw, and
prepare initial state vectors.

### Inputs

- `m::AbstractModel`: model object
- `df::DataFrame`: historical data. If `cond_type in [:semi, :full]`, then the
   final row of `df` should be the period containing conditional data

### Keyword Arguments

- `input_type::Symbol`: See documentation for `forecast_all`. Defaults to
  `:mode`
- `cond_type::Symbol`: See documentation for `forecast_all`. Defaults to `:none`
- `subset_inds::Vector{Int}`: indices specifying the draws we want to use. See
  `forecast_one` for more detail.
- `verbose::Symbol`: desired frequency of function progress messages printed to
  standard out. One of `:none`, `:low`, or `:high`. If `:low` or greater, prints
  location of input file.
- `procs::Vector{Int}`: list of worker processes that have been
  previously added by the user. Defaults to `[myid()]`

### Outputs

- `systems::DVector{System{Float64}}`: vector of `n_sim_forecast` many `System`
  objects, one for each draw
- `states::DVector{Vector{Float64}}`: vector of `n_sim` many final historical
  state vectors

### Notes

`prepare_forecast_inputs` calls `load_draws`, `prepare_systems`, and
  `prepare_states`. See those functions for thorough documentation.
"""
function prepare_forecast_inputs(m::AbstractModel, df::DataFrame;
    input_type::Symbol = :mode, cond_type::Symbol = :none,
    subset_inds::Vector{Int} = Vector{Int}(), verbose::Symbol = :low,
    procs::Vector{Int} = [myid()])

    # Reset procs to [myid()] if necessary
    procs = reset_procs(m, procs, Nullable(input_type))

    # Set up infiles
    params, TTT, RRR, CCC = load_draws(m, input_type; subset_inds = subset_inds,
                                      verbose = verbose, procs = procs)

    # Populate systems vector
    systems = prepare_systems(m, input_type, params, TTT, RRR, CCC; procs = procs)

    # Populate states vector
    kals = filter_all(m, df, systems; cond_type = cond_type, procs = procs)

    return systems, kals
end

"""
```
forecast_one(m, df; input_type = :mode, cond_type = :none, output_vars = [],
    subset_inds = [], verbose = :low, procs = [myid()])
```

Compute, save, and return `output_vars` for input draws given by `input_type`
and conditional data case given by `cond_type`.

### Inputs

- `m::AbstractModel`: model object
- `df::DataFrame`: Historical data. If `cond_type in [:semi, :full]`, then the
   final row of `df` should be the period containing conditional data.

### Keyword Arguments

- `input_type::Symbol`: See documentation for `forecast_all`. Defaults to
  `:mode`
- `cond_type::Symbol`: See documentation for `forecast_all`. Defaults to `:none`
- `output_vars::Vector{Symbol}`: vector of desired output variables. See Outputs
  section
- `subset_inds::Vector{Int}`: indices specifying the draws we want to use. If a
  more sophisticated selection criterion is desired, the user is responsible for
  determining the indices corresponding to that criterion. If `input_type` is
  not `subset`, `subset_inds` will be ignored.
- `subset_string::AbstractString`: short string identifying the subset to be
  appended to the output filenames. If `input_type = :subset` and
  `subset_string` is empty, an error is thrown.
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
  - `:shockdecstates`: `DArray{Float64, 4}` of state shock decompositions
  - `:shockdecobs`: `DArray{Float64, 4}` of observable shock decompositions
  - `:shockdecpseudo`: `DArray{Float64, 4}` of pseudo-observable shock
    decompositions (if a pseudo-measurement equation has been provided for this
    model type)
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
function forecast_one(m::AbstractModel{Float64}, df::DataFrame;
    input_type::Symbol = :mode, cond_type::Symbol = :none,
    output_vars::Vector{Symbol} = Vector{Symbol}(),
    subset_inds::Vector{Int} = Vector{Int}(),
    subset_string::AbstractString = "", verbose::Symbol = :low,
    procs::Vector{Int} = [myid()])

    ### 1. Setup

    # Reset procs to [myid()] if necessary
    procs = reset_procs(m, procs, Nullable(input_type))

    # Compute everything that will be needed to plot original output_vars
    output_vars = add_requisite_output_vars(output_vars)

    # Prepare forecast outputs
    forecast_output = Dict{Symbol, DArray{Float64}}()
    forecast_output_files = get_output_files(m, "forecast", input_type, output_vars, cond_type;
                                subset_string = subset_string)
    output_dir = rawpath(m, "forecast")

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println()
        info("Forecasting input_type = $input_type, cond_type = $cond_type...")
        println("Start time: $(now())")
        println("Forecast outputs will be saved in $output_dir")
    end

    # Prepare forecast inputs
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("\nPreparing forecast inputs...")
    end
    @time_verbose systems, kals = prepare_forecast_inputs(m, df; input_type = input_type,
            cond_type = cond_type, subset_inds = subset_inds, verbose = verbose,
            procs = procs)

    nprocs = length(procs)
    ndraws = length(systems)

    # Set forecast_pseudoobservables properly
    for output in output_vars
        if contains(string(output), "pseudo")
            m <= Setting(:forecast_pseudoobservables, true)
            break
        end
    end


    ### 2. Smoothed Histories

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
        @time_verbose histstates, histshocks, histpseudo =
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

        write_forecast_outputs(m, hists_to_compute, forecast_output_files, forecast_output; verbose = verbose)
    end


    ### 3. Forecasts

    forecast_vars = [:forecaststates, :forecastobs, :forecastpseudo, :forecastshocks]
    forecasts_to_compute = intersect(output_vars, forecast_vars)
    if !isempty(forecasts_to_compute)
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            println("\nForecasting $(forecasts_to_compute)...")
        end
        @time_verbose forecaststates, forecastobs, forecastpseudo, forecastshocks =
            forecast(m, systems, kals; cond_type = cond_type, procs = procs)

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

        write_forecast_outputs(m, forecasts_to_compute, forecast_output_files, forecast_output; verbose = verbose)
    end


    ### 4. Shock Decompositions
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

        write_forecast_outputs(m, shockdecs_to_compute, forecast_output_files, forecast_output; verbose = verbose)
    end


    ### 5. Trend

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

        write_forecast_outputs(m, trends_to_compute, forecast_output_files, forecast_output; verbose = verbose)
    end

    ### 6. Deterministic Trend

    dettrends_to_compute = intersect(output_vars, dettrend_vars)
    if !isempty(dettrends_to_compute)

        # Extract first smoothed historical state vectors for each
        # draw. This code will keep everything on its local node
        # rather than copying across nodes. The result is a
        # DVector{Vector{Float64}}.
        z0s = DArray(I->[convert(Vector{Float64}, slice(histstates, i, :, 1))::Vector{Float64} for i in first(I)],
                     (ndraws,), procs, [nprocs])

        # Compute deterministic trend
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            println("\nComputing deterministic trend for $(dettrends_to_compute)...")
        end
        @time_verbose dettrendstates, dettrendobs, dettrendpseudo =
            deterministic_trends(m, systems, z0s; procs = procs)

        forecast_output[:dettrendstates] = dettrendstates
        forecast_output[:dettrendobs]    = dettrendobs
        if :dettrendpseudo in output_vars
            forecast_output[:dettrendpseudo] = dettrendpseudo
        end

        dettrend_write_vars = [symbol("dettrend$c") for c in unique(map(get_class, dettrends_to_compute))]
        write_forecast_outputs(m, dettrends_to_compute, forecast_output_files, forecast_output; verbose = verbose)
    end


    # Return only saved elements of dict
    filter!((k, v) -> k in output_vars, forecast_output)
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("\nForecast complete: $(now())")
    end

    return forecast_output
end
