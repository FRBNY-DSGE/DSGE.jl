"""
```
forecast_all(m, cond_types, input_types, output_types)
```

Compute forecasts for all specified combinations of conditional data, input
types, and output types.

### Inputs

- `m::AbstractModel`: model object

- `cond_types::Vector{Symbol}`: conditional data types, any combination of

    - `:none`: no conditional data
    - `:semi`: use \"semiconditional data\" - average of quarter-to-date
      observations for high frequency series
    - `:full`: use \"conditional data\" - semiconditional plus nowcasts for
      desired observables

- `input_types::Vector{Symbol}`: which set of parameters to use, any combination of

    - `:mode`: forecast using the modal parameters only
    - `:mean`: forecast using the mean parameters only
    - `:init`: forecast using the initial parameter values only
    - `:full`: forecast using all parameters (full distribution)
    - `:subset`: forecast using a well-defined user-specified subset of draws

- `output_vars::Vector{Symbol}`: vector of desired output variables. See
  `forecast_one` for documentation of all possible `output_vars`

### Outputs

None. Output is saved to files returned by
`get_output_files(m, input_type, output_vars, cond_type)`
for each combination of `input_type` and `cond_type`.
"""
function forecast_all(m::AbstractModel,
                      cond_types::Vector{Symbol}  = Vector{Symbol}(),
                      input_types::Vector{Symbol} = Vector{Symbol}(),
                      output_vars::Vector{Symbol} = Vector{Symbol}();
                      procs::Vector{Int}          = [myid()])

    for cond_type in cond_types
        df = load_data(m; cond_type=cond_type, try_disk=true, verbose=:none)
        for input_type in input_types
            my_procs = if input_type in [:init, :mode, :mean]
                [myid()]
            else
                procs
            end
            forecast_one(m, df; cond_type=cond_type, input_type=input_type, output_vars=output_vars, procs=my_procs)
        end
    end
end

"""
`load_draws(m, input_type)`

Load and return parameter draws, transition matrices, and final state vectors
from Metropolis-Hastings. Single draws are reshaped to have additional singleton
dimensions, and missing variables without sufficient information are initialized
to null values of appropriate types.

### Inputs

- `m::AbstractModel`: model object
- `input_type::Symbol`: one of the options for `input_type` described in the
  documentation for `forecast_all`

### Outputs

- `params::Matrix{Float64}`: matrix of size `nsim` x `nparams` of parameter
  draws
- `TTT::Array{Float64, 3}`: array of size `nsim` x `nstates` x `nstates` of
  transition matrix draws
- `RRR::Array{Float64, 3}`: array of size `nsim` x `nstates` x `nshocks` of
  transition matrix draws
- `CCC::Array{Float64, 3}`: array of size `nsim` x `nstates` x 1 of transition
  matrix draws
- `zend::Matrix{Float64}`: matrix of size `nsim` x `nstates` of final historical
  period state vector draws

where `nsim` is the number of draws saved in Metropolis-Hastings.
"""
function load_draws(m::AbstractModel, input_type::Symbol)

    input_file_name = get_input_file(m, input_type)

    # Read infiles and set n_sim based on input_type type
    if input_type in [:mean, :mode]
        tmp = h5open(input_file_name, "r") do f
            map(Float64, read(f, "params"))
        end
        params = reshape(tmp, 1, size(tmp,1))
        TTT  = Array{Float64}(0,0,0)
        RRR  = Array{Float64}(0,0,0)
        CCC  = Array{Float64}(0,0,0)
        zend = Array{Float64}(0,0)
    elseif input_type in [:full]
        params, TTT, RRR, CCC, zend = h5open(input_file_name, "r") do f
            params = map(Float64, read(f, "mhparams"))
            TTT    = map(Float64, read(f, "mhTTT"))
            RRR    = map(Float64, read(f, "mhRRR"))
            zend   = map(Float64, read(f, "mhzend"))
            if "mhCCC" in names(f)
                CCC = map(Float64, read(f, "mhCCC"))
            else
                CCC = Array{Float64}(0,0,0)
            end
            params, TTT, RRR, CCC, zend
        end
    elseif input_type in [:init]
        init_parameters!(m)
        tmp = Float64[α.value for α in m.parameters]
        params = reshape(tmp, 1, size(tmp,1))
        TTT  = Array{Float64}(0,0,0)
        RRR  = Array{Float64}(0,0,0)
        CCC  = Array{Float64}(0,0,0)
        zend = Array{Float64}(0,0)
    end

    return params, TTT, RRR, CCC, zend
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
    elseif input_type in [:full]
        empty = isempty(CCC)
        nprocs = length(procs);
        systems = DArray((n_sim_forecast,), procs, [nprocs]) do I
            draw_inds = first(I)
            ndraws_local = Int(n_sim_forecast / nprocs)
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
prepare_states(m, input_type, cond_type, systems, params, df, zend;
    procs = [myid()])
```

Prepare the final historical state vector(s) for this combination of
forecast inputs. The final state vector is assumed to be from the period
corresponding to the final row of `df`. In the single-draw and
conditional cases, the final state vector is computed by running the
Kalman filter. In the multi-draw, unconditional cases, where the final
state vector has been successfully precomputed, final state vectors
are repackaged for inputs to the forecast.

### Inputs

- `m::AbstractModel`: model object
- `input_type::Symbol`: See documentation for `forecast_all`
- `cond_type::Symbol`: See documentation for `forecast_all`
- `systems::DVector{System{Float64}}`: vector of `n_sim_forecast` many `System`
  objects
- `params::Matrix{Float64}`: matrix of size `n_sim` x `nstates` of parameter
  draws from estimation. `systems[i]` is the `System` computed using
  `params[i*jstep, :]`
- `df::DataFrame`: historical data. If `cond_type in [:semi, :full]`, then the
   final row of `df` should be the period containing conditional data
- `zend::Matrix{Float64}`: matrix, either empty or of size `n_sim` x `nstates`,
  of final state vectors read in from `load_draws`

where `n_sim_forecast = n_sim / jstep` is the number of draws after thinning a
second time.

### Keyword arguments

- `procs::Vector{Int}`: list of worker processes that have been
  previously added by the user. Defaults to `[myid()]`

### Outputs

- `states::DVector{Vector{Float64}}`: vector of `n_sim` many final historical
  state vectors

### Notes

In all cases, the initial values in the forecast are the final historical state
vectors that are returned from the Kalman filter, not the Kalman or
Durbin-Koopman smoother (though these are the same in the case of the Kalman
smoother). These are the correct values to use because the Kalman filter returns
the actual mean (`s_{T|T}`) and variance (`P_{T|T}`) of the states given the
parameter draw, while the simulation smoother takes into account the uncertainty
in the parameter draw. Since we only compute one final historical state vector
for each parameter draw, we want to use the analytical mean and variance as the
starting points in the forecast.
"""
function prepare_states(m::AbstractModel, input_type::Symbol, cond_type::Symbol,
    systems::DVector{System{Float64}, Vector{System{Float64}}},
    params::Matrix{Float64}, df::DataFrame, zend::Matrix{Float64};
    procs::Vector{Int} = [myid()])

    # Reset procs to [myid()] if necessary
    procs = reset_procs(m, procs, Nullable(input_type))

    # Setup
    n_sim_forecast = length(systems)
    n_sim = size(params, 1)
    jstep = convert(Int, n_sim/n_sim_forecast)

    # If we just have one draw of parameters in mode, mean, or init case, then we don't have the
    # pre-computed system matrices. We now recompute them here by running the Kalman filter.
    if input_type in [:mean, :mode, :init]
        update!(m, vec(params))
        kal = filter(m, df, systems[1]; cond_type = cond_type, allout = true)
        # `kals` is a vector of length 1
        states = dfill(kal[:filt][:, end], (1,), procs);

    # If we have many draws, then we must package them into a vector of System objects.
    elseif input_type in [:full]
        if cond_type in [:none]
            nprocs = length(procs)
            states = DArray((n_sim_forecast,), procs, [nprocs]) do I
                draw_inds = first(I)
                ndraws_local = Int(n_sim_forecast / nprocs)
                localpart = Vector{Vector{Float64}}(ndraws_local)

                for i in draw_inds
                    j = i * jstep
                    i_local = mod(i-1, ndraws_local) + 1

                    localpart[i_local] = vec(zend[j, :])
                end
                return localpart
            end
        elseif cond_type in [:semi, :full]
            # We will need to re-run the entire filter/smoother so we can't do anything
            # here. The reason is that while we have $s_{T|T}$ we don't have $P_{T|T}$ and
            # thus can't "restart" the Kalman filter for the conditional data period.
            states = dfill(Vector{Float64}(), (0,), procs)
        else
            throw(ArgumentError("Not implemented for this `cond_type`."))
        end
    else
        throw(ArgumentError("Not implemented for this `input_type`."))
    end

    return states
end

"""
```
prepare_forecast_inputs(m, df; input_type = :mode, cond_type = :none,
    procs = [myid()])
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
- `procs::Vector{Int}`: list of worker processes that have been
  previously added by the user. Defaults to `[myid()]`

### Outputs

- `systems::DVector{System{Float64}}`: vector of `n_sim_forecast` many `System`
  objects, one for each draw
- `states::DVector{Vector{Float64}}`: vector of `n_sim` many final historical
  state vectors

### Notes

`prepare_forecast_inputs` calls `prepare_systems` and `prepare_states`. See
  those functions for thorough documentation.
"""
function prepare_forecast_inputs(m::AbstractModel, df::DataFrame;
    input_type::Symbol = :mode, cond_type::Symbol = :none,
    procs::Vector{Int} = [myid()])

    # Reset procs to [myid()] if necessary
    procs = reset_procs(m, procs, Nullable(input_type))

    # Set up infiles
    params, TTT, RRR, CCC, zend = load_draws(m, input_type)

    n_sim = size(params,1)
    jstep = get_jstep(m, n_sim)
    n_sim_forecast = convert(Int, n_sim/jstep)

    # Populate systems vector
    systems = prepare_systems(m, input_type, params, TTT, RRR, CCC; procs = procs)

    # Populate states vector
    states = prepare_states(m, input_type, cond_type, systems, params, df, zend; procs = procs)

    return systems, states
end

"""
```
forecast_one(m, df; input_type = :mode, cond_type = :none, output_vars = [],
    verbose = :low, procs = [myid()])
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
- `procs::Vector{Int}`: list of worker processes that have been previously added
  by the user. Defaults to `[myid()]`

### Output

- `forecast_outputs::Dict{Symbol, DArray{Float64}}`: dictionary of forecast
  outputs. Keys are `output_vars`, which is some subset of:

  + `:histstates::DArray{Float64, 3}`: smoothed historical states
  + `:histobs::DArray{Float64, 3}`: smoothed historical data
  + `:histpseudo::DArray{Float64, 3}`: smoothed historical pseudo-observables
    (if a pseudomeasurement equation has been provided for this model type)
  + `:histshocks::DArray{Float64, 3}`: smoothed historical shocks
  + `:forecaststates::DArray{Float64, 3}`: forecasted states
  + `:forecastobs::DArray{Float64, 3}`: forecasted observables
  + `:forecastpseudo::DArray{Float64, 3}`: forecasted pseudo-observables (if a
    pseudomeasurement equation has been provided for this model type)
  + `:forecastshocks::DArray{Float64, 3}`: forecasted shocks
  + `:shockdecstates::DArray{Float64, 4}`: state shock decompositions
  + `:shockdecobs::DArray{Float64, 4}`: observable shock decompositions
  + `:shockdecpseudo::DArray{Float64, 4}`: pseudo-observable shock
    decompositions (if a pseudomeasurement equation has been provided for this
    model type)

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
    output_vars::Vector{Symbol} = [], verbose::Symbol = :low,
    procs::Vector{Int} = [myid()])

    ### 1. Setup

    # Reset procs to [myid()] if necessary
    procs = reset_procs(m, procs, Nullable(input_type))

    # Prepare forecast inputs
    systems, states = prepare_forecast_inputs(m, df; input_type = input_type,
        cond_type = cond_type, procs = procs)

    nprocs = length(procs)
    ndraws = length(systems)

    # Prepare forecast outputs
    forecast_output = Dict{Symbol, DArray{Float64}}()
    forecast_output_files = get_output_files(m, input_type, output_vars, cond_type)
    output_dir = rawpath(m, "forecast")

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("\nForecasting input_type = $input_type, cond_type = $cond_type...")
        println("Forecast outputs will be saved in $output_dir")
    end

    # Inline definition s.t. the dicts forecast_output and forecast_output_files are accessible
    function write_forecast_outputs(vars::Vector{Symbol})
        for var in vars
            filepath = forecast_output_files[var]
            jldopen(filepath, "w") do file
                write_forecast_metadata(m, file, var)
            end
            write_darray(filepath, forecast_output[var])

            if VERBOSITY[verbose] >= VERBOSITY[:high]
                println(" * Wrote $(basename(filepath))")
            end
        end
    end


    ### 2. Smoothed Histories

    # Set forecast_pseudoobservables properly
    for output in output_vars
        if contains(string(output), "pseudo")
            m <= Setting(:forecast_pseudoobservables, true)
            break
        end
    end

    # Must re-run filter/smoother for conditional data in addition to explicit cases
    hist_vars = [:histstates, :histpseudo, :histshocks]
    shockdec_vars = [:shockdecstates, :shockdecpseudo, :shockdecobs]
    filterandsmooth_vars = vcat(hist_vars, shockdec_vars)

    if !isempty(intersect(output_vars, filterandsmooth_vars)) || cond_type in [:semi, :full]

        histstates, histshocks, histpseudo, zends =
            filterandsmooth_all(m, df, systems; cond_type = cond_type, procs = procs)

        # For conditional data, transplant the obs/state/pseudo vectors from hist to forecast
        if cond_type in [:semi, :full]
            T = n_mainsample_periods(m)

            forecast_output[:histstates] = convert(DArray, histstates[1:end, 1:end, 1:T])
            forecast_output[:histshocks] = convert(DArray, histshocks[1:end, 1:end, 1:T])
	        if :histpseudo in output_vars
                forecast_output[:histpseudo] = convert(DArray, histpseudo[1:end, 1:end, 1:T])
            end
        else
            forecast_output[:histstates] = histstates
            forecast_output[:histshocks] = histshocks
            if :histpseudo in output_vars
                forecast_output[:histpseudo] = histpseudo
            end
        end

        write_forecast_outputs(intersect(output_vars, hist_vars))
    end


    ### 3. Forecasts

    # For conditional data, use the end of the hist states as the initial state
    # vector for the forecast
    if cond_type in [:semi, :full]
        states = zends
    end

    forecast_vars = [:forecaststates, :forecastobs, :forecastpseudo, :forecastshocks]

    if !isempty(intersect(output_vars, forecast_vars))
        forecaststates, forecastobs, forecastpseudo, forecastshocks =
            forecast(m, systems, states; procs = procs)

        # For conditional data, transplant the obs/state/pseudo vectors from hist to forecast
        if cond_type in [:semi, :full]
            T = n_mainsample_periods(m)

            # Copy history of observables to make correct size
            histobs_cond = df_to_matrix(m, df; cond_type = cond_type)[:, index_mainsample_start(m)+T:end]
            histobs_cond = reshape(histobs_cond, (1, size(histobs_cond)...))

            nperiods = (size(histstates, 3) - T) + size(forecaststates, 3)
            nobs = n_observables(m)

            function cat_conditional(histvars::DArray{Float64, 3}, forecastvars::DArray{Float64, 3})
                nvars = size(histvars, 2)
                return DArray((ndraws, nvars, nperiods), procs, [nprocs, 1, 1]) do I
                    draw_inds = first(I)
                    hist_cond = convert(Array, histvars[draw_inds, :, T+1:end])
                    forecast  = convert(Array, forecastvars[draw_inds, :, :])
                    return cat(3, hist_cond, forecast)
                end
            end

            forecast_output[:forecaststates] = cat_conditional(histstates, forecaststates)
            forecast_output[:forecastshocks] = cat_conditional(histshocks, forecastshocks)
	        if :forecastpseudo in output_vars
                forecast_output[:forecastpseudo] = cat_conditional(histpseudo, forecastpseudo)
            end
            forecast_output[:forecastobs]    =
                DArray((ndraws, nobs, nperiods), procs, [nprocs, 1, 1]) do I
                    draw_inds = first(I)
                    ndraws_local = length(draw_inds)
                    hist_cond = repeat(histobs_cond, outer = [ndraws_local, 1, 1])
                    forecast  = convert(Array, forecastobs[draw_inds, :, :])
                    return cat(3, hist_cond, forecast)
                end
        else
            forecast_output[:forecaststates] = forecaststates
            forecast_output[:forecastshocks] = forecastshocks
            forecast_output[:forecastobs]    = forecastobs
            if :forecastpseudo in output_vars
                forecast_output[:forecastpseudo] = forecastpseudo
            end
        end

        write_forecast_outputs(intersect(output_vars, forecast_vars))
    end


    ### 4. Shock Decompositions

    if !isempty(intersect(output_vars, shockdec_vars))
        shockdecstates, shockdecobs, shockdecpseudo =
            shock_decompositions(m, systems, histshocks; procs = procs)

        forecast_output[:shockdecstates] = shockdecstates
        forecast_output[:shockdecobs]    = shockdecobs
        if :forecastpseudo in output_vars
            forecast_output[:shockdecpseudo] = shockdecpseudo
        end

        write_forecast_outputs(intersect(output_vars, shockdec_vars))
    end

    # Return only saved elements of dict
    filter!((k, v) -> k in output_vars, forecast_output)
    return forecast_output
end
