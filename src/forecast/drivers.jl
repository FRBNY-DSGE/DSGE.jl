"""
```
prepare_forecast_inputs!(m, input_type, cond_type, output_vars;
    df = DataFrame(), verbose = :none)
```

Add required outputs using `add_requisite_output_vars` and load data if
necessary.

### Inputs

- `m::AbstractModel`: model object
- `input_type::Symbol`: See `?forecast_one`.
- `cond_type::Symbol`: See `?forecast_one`.
- `output_vars::Vector{Symbol}`: vector of desired output variables. See
  `?forecast_one_draw`

### Keyword Arguments

- `df::DataFrame`: historical data. If `cond_type in [:semi, :full]`, then the
   final row of `df` should be the period containing conditional data. If not
   provided, then `df` will be loaded using `load_data` with the appropriate
   `cond_type`
- `verbose::Symbol`: desired frequency of function progress messages printed to
  standard out. One of `:none`, `:low`, or `:high`

### Outputs

- `output_vars`
- `df`
"""
function prepare_forecast_inputs!{S<:AbstractFloat}(m::AbstractModel{S},
    input_type::Symbol, cond_type::Symbol, output_vars::Vector{Symbol};
    df::DataFrame = DataFrame(), subset_inds::Range{Int64} = 1:0,
    verbose::Symbol = :none)

    # Compute everything that will be needed to plot original output_vars
    output_vars = add_requisite_output_vars(output_vars)

    # Get products and classes computed
    output_prods   = unique(map(get_product, output_vars))
    output_classes = unique(map(get_class,   output_vars))

    # Set forecast_pseudoobservables properly
    if any(class -> class == :pseudo, output_classes)
        m <= Setting(:forecast_pseudoobservables, true)
    end

    # Throw error if input_type = :subset but no subset_inds provided
    if input_type == :subset && isempty(subset_inds)
        error("Must supply nonempty subset_inds if input_type = :subset")
    end

    # Determine if we are only running IRFs. If so, we won't need to load data
    # below
    irfs_only = all(prod -> prod == :irf, output_prods)

    # Load data if not provided, else check data well-formed
    if !irfs_only
        if isempty(df)
            data_verbose = verbose == :none ? :none : :low
            df = load_data(m; cond_type = cond_type, try_disk = true, verbose = data_verbose)
        else
            @assert df[1, :date] == date_presample_start(m)
            @assert df[end, :date] == (cond_type == :none ? date_mainsample_end(m) : date_conditional_end(m))
        end
    end

    return output_vars, df
end

"""
```
load_draws(m, input_type; subset_inds = 1:0, verbose = :low)

load_draws(m, input_type, block_inds; verbose = :low)
```

Load and return parameter draws from Metropolis-Hastings.

### Inputs

- `m::AbstractModel`: model object
- `input_type::Symbol`: one of the options for `input_type` described in the
  documentation for `forecast_one`
- `block_inds::Range{Int64}`: indices of the current block (already indexed by
  `jstep`) to be read in. Only used in second method

### Keyword Arguments

- `subset_inds::Range{Int64}`: indices specifying the subset of draws to be read
  in. Only used in first method
- `verbose::Symbol`: desired frequency of function progress messages printed to
  standard out. One of `:none`, `:low`, or `:high`. If `:low` or greater, prints
  location of input file.

### Outputs

- `params`: first method returns a single parameter draw of type
  `Vector{Float64}`. Second method returns a `Vector{Vector{Float64}}` of
  parameter draws for this block.
"""
function load_draws(m::AbstractModel, input_type::Symbol; subset_inds::Range{Int64} = 1:0,
    verbose::Symbol = :low)

    input_file_name = get_forecast_input_file(m, input_type)
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("Loading draws from $input_file_name")
    end

    # Load single draw
    if input_type in [:mean, :mode]

        params = convert(Vector{Float64}, h5read(input_file_name, "params"))

    # Load full distribution
    elseif input_type == :full

        params = map(Float64, h5read(input_file_name, "mhparams"))

    # Load subset of full distribution
    elseif input_type == :subset

        if isempty(subset_inds)
            error("Must supply nonempty range of subset_inds if input_type == :subset")
        else
            params = map(Float64, h5read(input_file_name, "mhparams", (subset_inds, :)))
        end

    # Return initial parameters of model object
    elseif input_type == :init

        init_parameters!(m)
        tmp = map(α -> α.value, m.parameters)
        params = convert(Vector{Float64}, tmp)

    end

    return params
end

function load_draws(m::AbstractModel, input_type::Symbol, block_inds::Range{Int64};
                    verbose::Symbol = :low)

    input_file_name = get_forecast_input_file(m, input_type)
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("Loading draws from $input_file_name")
    end

    if input_type in [:full, :subset]
        if isempty(block_inds)
            error("Must supply nonempty range of block_inds for this load_draws method")
        else
            ndraws = length(block_inds)
            params = Vector{Vector{Float64}}(ndraws)
            for (i, j) in zip(1:ndraws, block_inds)
                params[i] = vec(map(Float64, h5read(input_file_name, "mhparams", (j, :))))
            end
            return params
        end
    else
        error("This load_draws method can only be called with input_type in [:full, :subset]")
    end

end

"""
```
forecast_one(m, input_type, cond_type, output_vars; df = DataFrame(),
    subset_inds = 1:0, forecast_string = "", verbose = :low)
```

Compute, save, and return `output_vars` for input draws given by `input_type`
and conditional data case given by `cond_type`.

### Inputs

- `m::AbstractModel`: model object

- `input_type::Symbol`: one of:

```
  - `:mode`: forecast using the modal parameters only
  - `:mean`: forecast using the mean parameters only
  - `:init`: forecast using the initial parameter values only
  - `:full`: forecast using all parameters (full distribution)
  - `:subset`: forecast using a well-defined user-specified subset of draws
```

- `cond_type::Symbol`: one of:

```
  - `:none`: no conditional data
  - `:semi`: use \"semiconditional data\" - average of quarter-to-date
    observations for high frequency series
  - `:full`: use \"conditional data\" - semiconditional plus nowcasts for
    desired observables
```

- `output_vars::Vector{Symbol}`: vector of desired output variables. See
  `?forecast_one_draw`.

### Keyword Arguments

- `df::DataFrame`: Historical data. If `cond_type in [:semi, :full]`, then the
   final row of `df` should be the period containing conditional data. If not
   provided, will be loaded using `load_data` with the appropriate `cond_type`
- `subset_inds::Range{Int64}`: indices specifying the draws we want to use. If a
  more sophisticated selection criterion is desired, the user is responsible for
  determining the indices corresponding to that criterion. If `input_type` is
  not `subset`, `subset_inds` will be ignored
- `forecast_string::AbstractString`: short string identifying the subset to be
  appended to the output filenames. If `input_type = :subset` and
  `forecast_string` is empty, an error is thrown.
- `verbose::Symbol`: desired frequency of function progress messages printed to
  standard out. One of `:none`, `:low`, or `:high`.

### Outputs

None. Output is saved to files returned by
`get_forecast_output_files(m, input_type, cond_type, output_vars)`.
"""
function forecast_one(m::AbstractModel{Float64},
    input_type::Symbol, cond_type::Symbol, output_vars::Vector{Symbol};
    df::DataFrame = DataFrame(), subset_inds::Range{Int64} = 1:0,
    forecast_string::AbstractString = "", verbose::Symbol = :low)

    ### Common Setup

    # Add necessary output_vars and load data
    output_vars, df = prepare_forecast_inputs!(m, input_type, cond_type, output_vars;
                                               df = df, verbose = verbose)

    # Get output file names
    forecast_output = Dict{Symbol, Array{Float64}}()
    forecast_output_files = get_forecast_output_files(m, input_type, cond_type, output_vars;
                                                      forecast_string = forecast_string)
    output_dir = rawpath(m, "forecast")

    # Print
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        info("Forecasting input_type = $input_type, cond_type = $cond_type...")
        println("Start time: $(now())")
        println("Forecast outputs will be saved in $output_dir")
    end


    ### Single-Draw Forecasts

    if input_type in [:mode, :mean, :init]

        tic()

        params = load_draws(m, input_type; verbose = verbose)
        forecast_output = forecast_one_draw(m, input_type, cond_type, output_vars,
                                            params, df, verbose = verbose)

        write_forecast_outputs(m, input_type, output_vars, forecast_output_files,
                               forecast_output; block_number = Nullable{Int64}(),
                               verbose = verbose)

        if VERBOSITY[verbose] >= VERBOSITY[:low]
            total_forecast_time     = toq()
            total_forecast_time_min = total_forecast_time/60

            println("\nTotal time to forecast: $total_forecast_time_min minutes")
        end


    ### Multiple-Draw Forecasts

    elseif input_type in [:full, :subset]

        # Block info
        block_inds = forecast_block_inds(m, input_type; subset_inds = subset_inds)
        nblocks = length(block_inds)
        start_block = isnull(forecast_start_block(m)) ? 1 : get(forecast_start_block(m))

        # Info needed for printing progress
        total_forecast_time = 0.0
        block_verbose = verbose == :none ? :none : :low

        for block = start_block:nblocks
            if VERBOSITY[verbose] >= VERBOSITY[:low]
                println()
                info("Forecasting block $block of $nblocks...")
            end
            tic()

            # Get to work!
            params = load_draws(m, input_type, block_inds[block]; verbose = verbose)

            forecast_outputs = pmap(param -> forecast_one_draw(m, input_type, cond_type, output_vars,
                                                               param, df, verbose = verbose),
                                    params)

            # If some element of forecast_outputs is a RemoteException, rethrow the exception
            ind_ex = findfirst(x -> isa(x, RemoteException), forecast_outputs)
            if ind_ex > 0
                ex = forecast_outputs[ind_ex].captured
                throw(ex)
            else
                forecast_outputs = convert(Vector{Dict{Symbol, Array{Float64}}}, forecast_outputs)
            end

            # Assemble outputs from this block and write to file
            forecast_output = assemble_block_outputs(forecast_outputs)
            write_forecast_outputs(m, input_type, output_vars, forecast_output_files,
                                   forecast_output; block_number = Nullable(block),
                                   verbose = block_verbose, block_inds = block_inds[block])
            gc()

            # Calculate time to complete this block, average block time, and
            # expected time to completion
            if VERBOSITY[verbose] >= VERBOSITY[:low]
                block_time = toq()
                total_forecast_time += block_time
                total_forecast_time_min     = total_forecast_time/60
                expected_time_remaining     = (total_forecast_time/block)*(nblocks - block)
                expected_time_remaining_min = expected_time_remaining/60

                println("\nCompleted $block of $nblocks blocks.")
                println("Total time to compute $block blocks: $total_forecast_time_min minutes")
                println("Expected time remaining in forecast: $expected_time_remaining_min minutes")
            end
        end # of loop through blocks

    end # of input_type

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("\nForecast complete: $(now())")
    end
end

"""
```
forecast_one_draw(m, input_type, cond_type, output_vars; params, df;
    verbose = :low)
```

Compute `output_vars` for a single parameter draw, `params`. Called by
`forecast_one`.

### Inputs

- `m::AbstractModel{Float64}`: model object
- `input_type::Symbol`: See `?forecast_one`.
- `cond_type::Symbol`: See `?forecast_one`.
- `output_vars::Vector{Symbol}`: vector of desired output variables. See Outputs
  section
- `params::Vector{Float64}`: parameter vector
- `df::DataFrame`: historical data.
- `verbose::Symbol`: desired frequency of function progress messages printed to
  standard out. One of `:none`, `:low`, or `:high`.

### Output

- `forecast_outputs::Dict{Symbol, Array{Float64}}`: dictionary of forecast
  outputs. Keys are `output_vars`, which is some subset of:

```
  - `:histstates`: `Matrix{Float64}` of smoothed historical states
  - `:histobs`: `Matrix{Float64}` of smoothed historical data
es  - `:histpseudo`: `Matrix{Float64}` of smoothed historical
    pseudo-observables (if a pseudo-measurement equation has been provided for
    this model type)
  - `:histshocks`: `Matrix{Float64}` of smoothed historical shocks
  - `:forecaststates`: `Matrix{Float64}` of forecasted states
  - `:forecastobs`: `Matrix{Float64}` of forecasted observables
  - `:forecastpseudo`: `Matrix{Float64}` of forecasted pseudo-observables
  - `:forecastshocks`: `Matrix{Float64}` of forecasted shocks
  - `:bddforecaststates`, `:bddforecastobs`, `:bddforecastpseudo`, and
    `:bddforecastshocks`: `Matrix{Float64}`s of forecasts where we enforce
    the zero lower bound to be `forecast_zlb_value(m)`
  - `:shockdecstates`: `Array{Float64, 3}` of state shock decompositions
  - `:shockdecobs`: `Array{Float64, 3}` of observable shock decompositions
  - `:shockdecpseudo`: `Array{Float64, 3}` of pseudo-observable shock
    decompositions
  - `:dettrendstates`: `Matrix{Float64}` of state deterministic trends
  - `:dettrendobs`: `Matrix{Float64}` of observable deterministic trends
  - `:dettrendpseudo`: `Matrix{Float64}` of pseudo-observable deterministic
    trends
  - `:trendstates`: `Vector{Float64}` of state trends, i.e. the `CCC` vector
  - `:trendobs`: `Vector{Float64}` of observable trends, i.e. the `DD` vector
  - `:trendpseudo`: `Vector{Float64}` of pseudo-observable trends, i.e. the
    `DD_pseudo` vector
  - `:irfstates`: `Array{Float64, 3}` of state impulse responses
  - `:irfobs`: `Array{Float64, 3}` of observable impulse responses
  - `:irfpseudo`: `Array{Float64, 3}` of pseudo-observable impulse responses
```
"""
function forecast_one_draw(m::AbstractModel{Float64}, input_type::Symbol, cond_type::Symbol,
    output_vars::Vector{Symbol}, params::Vector{Float64}, df::DataFrame; verbose::Symbol = :low)

    ### Setup

    # Are we only running IRFs?
    output_prods = map(get_product, output_vars)
    irfs_only = all(x -> x == :irf, output_prods)

    # Compute state space and run Kalman filter
    update!(m, params)
    system = compute_system(m)
    if !irfs_only
        kal = filter(m, df, system; cond_type = cond_type)
    end

    # Initialize dictionary
    forecast_output = Dict{Symbol, Array{Float64}}()


    ### 1. Smoothed Histories

    # Must run smoother for conditional data in addition to explicit cases
    hist_vars = [:histstates, :histpseudo, :histshocks]
    shockdec_vars = [:shockdecstates, :shockdecpseudo, :shockdecobs]
    dettrend_vars = [:dettrendstates, :dettrendpseudo, :dettrendobs]
    smooth_vars = vcat(hist_vars, shockdec_vars, dettrend_vars)
    hists_to_compute = intersect(output_vars, hist_vars)

    if !isempty(intersect(output_vars, smooth_vars)) || cond_type in [:semi, :full]

        histstates, histshocks, histpseudo, initial_states =
            smooth(m, df, system, kal; cond_type = cond_type)

        # For conditional data, transplant the obs/state/pseudo vectors from hist to forecast
        if cond_type in [:full, :semi]
            T = n_mainsample_periods(m)

            forecast_output[:histstates] = transplant_history(histstates, T)
            forecast_output[:histshocks] = transplant_history(histshocks, T)
            forecast_output[:histpseudo] = transplant_history(histpseudo, T)
        else
            forecast_output[:histstates] = histstates
            forecast_output[:histshocks] = histshocks
            forecast_output[:histpseudo] = histpseudo
        end
    end


    ### 2. Forecasts

    # 2A. Unbounded forecasts

    forecast_vars = [:forecaststates, :forecastobs, :forecastpseudo, :forecastshocks]
    forecasts_to_compute = intersect(output_vars, forecast_vars)

    if !isempty(forecasts_to_compute)
        forecaststates, forecastobs, forecastpseudo, forecastshocks =
            forecast(m, system, kal; cond_type = cond_type, enforce_zlb = false)

        # For conditional data, transplant the obs/state/pseudo vectors from hist to forecast
        if cond_type in [:full, :semi]
            forecast_output[:forecaststates] = transplant_forecast(histstates, forecaststates, T)
            forecast_output[:forecastshocks] = transplant_forecast(histshocks, forecastshocks, T)
            forecast_output[:forecastpseudo] = transplant_forecast(histpseudo, forecastpseudo, T)
            forecast_output[:forecastobs]    = transplant_forecast_observables(histstates, forecastobs, system, T)
        else
            forecast_output[:forecaststates] = forecaststates
            forecast_output[:forecastshocks] = forecastshocks
            forecast_output[:forecastpseudo] = forecastpseudo
            forecast_output[:forecastobs]    = forecastobs
        end
    end


    # 2B. Bounded forecasts

    forecast_vars_bdd = [:bddforecaststates, :bddforecastobs, :bddforecastpseudo, :bddforecastshocks]
    forecasts_to_compute = intersect(output_vars, forecast_vars_bdd)

    if !isempty(forecasts_to_compute)
        forecaststates, forecastobs, forecastpseudo, forecastshocks =
            forecast(m, system, kal; cond_type = cond_type, enforce_zlb = true)

        # For conditional data, transplant the obs/state/pseudo vectors from hist to forecast
        if cond_type in [:full, :semi]
            forecast_output[:bddforecaststates] = transplant_forecast(histstates, forecaststates, T)
            forecast_output[:bddforecastshocks] = transplant_forecast(histshocks, forecastshocks, T)
            forecast_output[:bddforecastpseudo] = transplant_forecast(histpseudo, forecastpseudo, T)
            forecast_output[:bddforecastobs]    = transplant_forecast_observables(histstates, forecastobs, system, T)
        else
            forecast_output[:bddforecaststates] = forecaststates
            forecast_output[:bddforecastshocks] = forecastshocks
            forecast_output[:bddforecastpseudo] = forecastpseudo
            forecast_output[:bddforecastobs]    = forecastobs
        end
    end


    ### 3. Shock Decompositions

    shockdecs_to_compute = intersect(output_vars, shockdec_vars)

    if !isempty(shockdecs_to_compute)
        shockdecstates, shockdecobs, shockdecpseudo = shock_decompositions(m, system, histshocks)

        forecast_output[:shockdecstates] = shockdecstates
        forecast_output[:shockdecobs]    = shockdecobs
        forecast_output[:shockdecpseudo] = shockdecpseudo
    end


    ### 4. Trend

    trend_vars = [:trendstates, :trendobs, :trendpseudo]
    trends_to_compute = intersect(output_vars, trend_vars)

    if !isempty(trends_to_compute)
        trendstates, trendobs, trendpseudo = trends(system)

        forecast_output[:trendstates] = trendstates
        forecast_output[:trendobs]    = trendobs
        forecast_output[:trendpseudo] = trendpseudo
    end


    ### 5. Deterministic Trend

    dettrends_to_compute = intersect(output_vars, dettrend_vars)

    if !isempty(dettrends_to_compute)
        dettrendstates, dettrendobs, dettrendpseudo = deterministic_trends(m, system, initial_states)

        forecast_output[:dettrendstates] = dettrendstates
        forecast_output[:dettrendobs]    = dettrendobs
        forecast_output[:dettrendpseudo] = dettrendpseudo
    end


    ### 6. Impulse Responses

    irf_vars = [:irfstates, :irfobs, :irfpseudo]
    irfs_to_compute = intersect(output_vars, irf_vars)

    if !isempty(irfs_to_compute)
        irfstates, irfobs, irfpseudo = impulse_responses(m, system)

        forecast_output[:irfstates] = irfstates
        forecast_output[:irfobs] = irfobs
        forecast_output[:irfpseudo] = irfpseudo
    end


    ### Return only desired output_vars

    for key in keys(forecast_output)
        if !(key in output_vars)
            delete!(forecast_output, key)
        end
    end
    return forecast_output
end