"""
```
forecast_all(m::AbstractModel, cond_types::Vector{Symbol}
input_types::Vector{Symbol}, output_types::Vector{Symbol}
```

Compute forecasts for all specified combinations of conditional data, input types, and
output types.

# Arguments

- `m`: model object
- `cond_types`: conditional data type, any combination of
    - `:none`: no conditional data
    - `:semi`: use \"semiconditional data\" - average of quarter-to-date
      observations for high frequency series
    - `:full`: use \"conditional data\" - semiconditional plus nowcasts for
      desired observables
- `input_types`: which set of parameters to use, any combination of
    - `:mode`: forecast using the modal parameters only
    - `:mean`: forecast using the mean parameters only
    - `:init`: forecast using the initial parameter values only
    - `:full`: forecast using all parameters (full distribution)
    - `:subset`: forecast using a well-defined user-specified subset of draws
- `output_types`: forecast routine outputs to compute, any combination of
    - `:states`: smoothed states (history) for all specified conditional data types
    - `:shocks`: smoothed shocks (history, standardized) for all specified
      conditional data types
    - `:shocks_nonstandardized`: smoothed shocks (history, non-standardized) for
      all specified conditional data types
    - `:forecast`: forecast of states and observables for all specified
      conditional data types, as well as shocks that produced them
    - `:shockdec`: shock decompositions (history) of states and observables for
      all specified conditional data types
    - `:dettrend`: deterministic trend (history) of states and observables for
      all specified conditional data types
    - `:counter`: counterfactuals (history) of states and observables for all
      specified conditional data types
    - `:simple`: smoothed states, forecast of states, forecast of observables
      for *unconditional* data only
    - `:all`: smoothed states (history), smoothed shocks (history, standardized), smoothed
      shocks (history, non-standardized), shock decompositions (history), deterministic
      trend (history), counterfactuals (history), forecast, forecast shocks drawn, shock
      decompositions (forecast), deterministic trend (forecast), counterfactuals (forecast)
   Note that some similar outputs may or may not fall under the \"forecast_all\" framework,
   including
    - `:mats`: recompute system matrices (TTT, RRR, CCC) given parameters only
    - `:zend`: recompute final state vector (s_{T}) given parameters only
    - `:irfs`: impulse response functions

Outputs
-------

- todo
"""
function forecast_all(m::AbstractModel,
                      cond_types::Vector{Symbol}   = Vector{Symbol}(),
                      input_types::Vector{Symbol}  = Vector{Symbol}(),
                      output_types::Vector{Symbol} = Vector{Symbol}())

    for cond_type in cond_types
        df = load_data(m; cond_type=cond_type, try_disk=true, verbose=:none)
        for input_type in input_types
            # Take the union of all output variables specified by output_types
            all_output_vars = map(x -> get_output_vars(m, x), output_types)
            output_vars = union(all_output_vars...)

            forecast_one(m, df; cond_type=cond_type, input_type=input_type, output_vars = output_vars)
        end
    end

end

"""
`load_draws(m, input_type)`

Load and return draws from Metropolis-Hastings, after some slight transformations. Single
draws are reshaped to have additional singleton dimensions, and missing variables without
sufficient information are initialized to null values of appropriate types.

### Outputs
- `params`: Matrix{Float64} of size (nsim, nparams)
- `TTT`: Array{Float64,3} of size (nsim, nequations, nstates)
- `RRR`: Array{Float64,3} of size (nsim, nequations, nshocks)
- `CCC`: Array{Float64,3} of size (nsim, nequations, 1)
- `zend`: Matrix{Float64} of size (nsim, nstates)
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
prepare_systems(m::AbstractModel, input_type::Symbol, params::Matrix{Float64},
TTT::Array{Float64,3}, RRR::Array{Float64,3}, CCC::Array{Float64,3})
```

Return Vector of System objects constructed from the given sampling outputs. In the one-draw
case (mode, mean, init), we recompute the entire system. In the many-draw case (full, or subset),
we package the outputs only. Recomputing the entire system in the many-draw case remains to
be implemented.
"""
function prepare_systems(m::AbstractModel, input_type::Symbol,
    params::Matrix{Float64}, TTT::Array{Float64, 3}, RRR::Array{Float64, 3},
    CCC::Array{Float64, 3}; my_procs::Vector{Int} = [myid()])

    # Setup
    n_sim = size(params,1)
    jstep = get_jstep(m, n_sim)
    n_sim_forecast = convert(Int, n_sim/jstep)

    if input_type in [:mean, :mode, :init]
        update!(m, vec(params))
        systems = dfill(compute_system(m), (1,), my_procs)
    elseif input_type in [:full]
        empty = isempty(CCC)
        nprocs = length(my_procs);
        systems = DArray((n_sim_forecast,), my_procs, [nprocs]) do I
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

                # Prepare system
                localpart[i_local] = System(trans_j, meas_j)
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
prepare_states(m::AbstractModel, input_type::Symbol, cond_type::Symbol,
               systems::Vector{System{Float64}}, params::Matrix{Float64}, df::DataFrame,
               zend::Matrix{Float64})
```

Return the final state vector(s) for this combination of inputs. The final state vector is
determined to be that s_{T} such that `T == size(df,1)`. Often, the final state vector is
computed by applying to the Kalman filter. In cases where the final state vector appears
to be successfully precomputed (such as full distribution input) but the data are
conditional data, then the final state vector is adjusted accordingly.

"""
function prepare_states(m::AbstractModel, input_type::Symbol, cond_type::Symbol,
    systems::DArray{System{Float64}, 1, Vector{System{Float64}}},
    params::Matrix{Float64}, df::DataFrame, zend::Matrix{Float64};
    my_procs::Vector{Int} = [myid()])

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
        states = dfill(kal[:filt][:, end], (1,), my_procs);

    # If we have many draws, then we must package them into a vector of System objects.
    elseif input_type in [:full]
        if cond_type in [:none]
            nprocs = length(my_procs)
            states = DArray((n_sim_forecast,), my_procs, [nprocs]) do I
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
            states = dfill(Vector{Float64}(), (0,), my_procs)
        end
    else
        throw(ArgumentError("Not implemented."))
    end

    return states
end

"""
```
prepare_forecast_inputs(m::AbstractModel, df::DataFrame; input_type::Symbol =
    :mode, cond_type::Symbol = :none)
```

Load draws for this input type, prepare a System object for each draw, and prepare initial
state vectors.
"""
function prepare_forecast_inputs(m::AbstractModel, df::DataFrame;
    input_type::Symbol = :mode, cond_type::Symbol = :none,
    my_procs::Vector{Int} = [myid()])

    # Set up infiles
    params, TTT, RRR, CCC, zend = load_draws(m, input_type)

    n_sim = size(params,1)
    jstep = get_jstep(m, n_sim)
    n_sim_forecast = convert(Int, n_sim/jstep)

    # Populate systems vector
    systems = prepare_systems(m, input_type, params, TTT, RRR, CCC; my_procs = my_procs)

    # Populate states vector
    states = prepare_states(m, input_type, cond_type, systems, params, df, zend; my_procs = my_procs)

    return systems, states
end

"""
```
forecast_one(m::AbstractModel, df::DataFrame; input_type::Symbol  = :mode,
    output_type::Symbol = :simple, cond_type::Symbol = :none)
```

Compute, save, and return forecast outputs given by `output_type` for input draws given by
`input_type` and conditional data case given by `cond_type`.

"""
function forecast_one(m::AbstractModel{Float64}, df::DataFrame;
    input_type::Symbol = :mode, cond_type::Symbol = :none,
    output_vars::Vector{Symbol} = [], verbose::Symbol = :low,
    my_procs::Vector{Int} = [myid()])

    ### 1. Setup

    # Use only one process for a single draw
    if input_type in [:init, :mode, :mean]
        my_procs = [myid()]
    end

    # Prepare forecast inputs
    systems, states = prepare_forecast_inputs(m, df; input_type = input_type,
        cond_type = cond_type, my_procs = my_procs)

    nprocs = length(my_procs)
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
            file = forecast_output_files[var]
            write_darray(file, forecast_output[var])
            if VERBOSITY[verbose] >= VERBOSITY[:high]
                println(" * Wrote $(basename(file))")
            end
        end
    end


    ### 2. Smoothed Histories

    # Set forecast_pseudoobservables properly
    for output in output_vars
        if contains(string(output), "pseudo")
            update!(m.settings[:forecast_pseudoobservables], Setting(:forecast_pseudoobservables, true))
            break
        end
    end

    # Must re-run filter/smoother for conditional data in addition to explicit cases
    hist_vars = [:histstates, :histpseudo, :histshocks]
    shockdec_vars = [:shockdecstates, :shockdecpseudo, :shockdecobs]
    filterandsmooth_vars = vcat(hist_vars, shockdec_vars)

    if !isempty(intersect(output_vars, filterandsmooth_vars)) || cond_type in [:semi, :full]

        histstates, histshocks, histpseudo, zends =
            filterandsmooth(m, df, systems; cond_type = cond_type, my_procs = my_procs)

        # For conditional data, transplant the obs/state/pseudo vectors from hist to forecast
        if cond_type in [:semi, :full]
            T = DSGE.subtract_quarters(date_forecast_start(m), date_prezlb_start(m))

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

        write_forecast_outputs(hist_vars)
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
            forecast(m, systems, states; my_procs = my_procs)

        # For conditional data, transplant the obs/state/pseudo vectors from hist to forecast
        if cond_type in [:semi, :full]
            T = DSGE.subtract_quarters(date_forecast_start(m), date_prezlb_start(m))

            # Copy history of observables to make correct size
            histobs_cond = df_to_matrix(m, df; cond_type = cond_type)[:, index_prezlb_start(m)+T:end]
            histobs_cond = reshape(histobs_cond, (1, size(histobs_cond)...))

            nperiods = (size(histstates, 3) - T) + size(forecaststates, 3)
            nobs = n_observables(m)

            function cat_conditional(histvars::DArray{Float64, 3}, forecastvars::DArray{Float64, 3})
                nvars = size(histvars, 2)
                return DArray((ndraws, nvars, nperiods), my_procs, [nprocs, 1, 1]) do I
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
                DArray((ndraws, nobs, nperiods), my_procs, [nprocs, 1, 1]) do I
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

        write_forecast_outputs(forecast_vars)
    end


    ### 4. Shock Decompositions

    if !isempty(intersect(output_vars, shockdec_vars))
        shockdecstates, shockdecobs, shockdecpseudo =
            shock_decompositions(m, systems, histshocks; my_procs = my_procs)

        forecast_output[:shockdecstates] = shockdecstates
        forecast_output[:shockdecobs]    = shockdecobs
        if :forecastpseudo in output_vars
            forecast_output[:shockdecpseudo] = shockdecpseudo
        end

        write_forecast_outputs(shockdec_vars)
    end

    # Return only saved elements of dict
    filter!((k, v) -> k in output_vars, forecast_output)
    return forecast_output
end
