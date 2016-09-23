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

function get_jstep(m, n_sim)
    if n_sim == 1
        jstep = 1
    else
        jstep = get_setting(m, :forecast_jstep)
    end
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
    systems::Vector{System{Float64}}, params::Matrix{Float64}, df::DataFrame,
    zend::Matrix{Float64})

    # Setup and preallocate
    n_sim_forecast = length(systems)
    n_sim = size(params,1)
    jstep = convert(Int, n_sim/n_sim_forecast)
    states = Vector{Vector{Float64}}(n_sim_forecast)

    # If we just have one draw of parameters in mode, mean, or init case, then we don't have the
    # pre-computed system matrices. We now recompute them here by running the Kalman filter.
    if input_type in [:mean, :mode, :init]
        update!(m, vec(params))
        kals = filter(m, df, systems; cond_type = cond_type, allout = true)
        # `kals` is a vector of length 1
        states[1] = kals[1][:filt][:, end]

    # If we have many draws, then we must package them into a vector of System objects.
    elseif input_type in [:full]
        if cond_type in [:none]
            # TODO if zend is empty for some reason, we should be able to recompute here
            for i in 1:n_sim_forecast
                j = i * jstep
                states[i] = vec(zend[j,:])
            end
        elseif cond_type in [:semi, :full]
            # We will need to re-run the entire filter/smoother so we can't do anything
            # here. The reason is that while we have $s_{T|T}$ we don't have $P_{T|T}$ and
            # thus can't "restart" the Kalman filter for the conditional data period.
            nothing
        end
    else
        throw(ArgumentError("Not implemented."))
    end

    return states
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
    params::Matrix{Float64}, TTT::Array{Float64,3}, RRR::Array{Float64,3},
    CCC::Array{Float64,3})

    # Setup and preallocate
    n_sim = size(params,1)
    jstep = get_jstep(m, n_sim)
    n_sim_forecast = convert(Int, n_sim/jstep)
    systems = Vector{System{Float64}}(n_sim_forecast)

    if input_type in [:mean, :mode, :init]
        update!(m, vec(params))
        systems[1] = compute_system(m)
    elseif input_type in [:full]
        empty = isempty(CCC)
        # TODO parallelize
        for i in 1:n_sim_forecast
            j = i * jstep
            # Prepare transition eq
            TTT_j  = squeeze(TTT[j,:,:],1)
            RRR_j  = squeeze(RRR[j,:,:],1)

            if empty
                trans_j = Transition(TTT_j, RRR_j)
            else
                CCC_j = squeeze(CCC[j,:,:],1)
                trans_j = Transition(TTT_j, RRR_j, CCC_j)
            end

            # Prepare measurement eq
            params_j = vec(params[j,:])
            update!(m, params_j)
            meas_j   = measurement(m, trans_j; shocks = true)

            # Prepare system
            systems[i] = System(trans_j, meas_j)
        end
    else
        throw(ArgumentError("Not implemented."))
    end

    return systems
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
                      input_type::Symbol  = :mode,
                      cond_type::Symbol   = :none)
    # Some variables
    n_states = n_states_augmented(m)

    # Set up infiles
    params, TTT, RRR, CCC, zend = load_draws(m, input_type)

    n_sim = size(params,1)
    jstep = get_jstep(m, n_sim)
    n_sim_forecast = convert(Int, n_sim/jstep)

    # Populate systems vector
    systems = prepare_systems(m, input_type, params, TTT, RRR, CCC)

    # Populate states vector
    states = prepare_states(m, input_type, cond_type, systems, params, df, zend)

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
function forecast_one(m::AbstractModel, df::DataFrame;
                      input_type::Symbol = :mode,
                      cond_type::Symbol = :none,
                      output_vars::Vector{Symbol} = [],
                      verbose::Symbol = :low)
    ### 1. Setup

    # Prepare forecast inputs
    systems, states = prepare_forecast_inputs(m, df; input_type = input_type,
        cond_type = cond_type)
    ndraws = length(systems)

    # Prepare forecast outputs
    forecast_output = Dict{Symbol, Vector{Array{Float64}}}()
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
            jldopen(file, "w") do f
                write(f, string(var), forecast_output[var])
                if VERBOSITY[verbose] >= VERBOSITY[:high]
                    println(" * Wrote $(basename(file))")
                end
            end
        end
    end


    ### 2. Smoothed Histories

    # must re-run filter/smoother for conditional data in addition to explicit cases
    hist_vars = [:histstates, :histpseudo, :histshocks]
    shockdec_vars = [:shockdecstates, :shockdecpseudo, :shockdecobs]
    filterandsmooth_vars = vcat(hist_vars, shockdec_vars)

    if !isempty(intersect(output_vars, filterandsmooth_vars)) || cond_type in [:semi, :full]

        histstates, histshocks, histpseudo, zends = filterandsmooth(m, df, systems; cond_type = cond_type)

        # For conditional data, transplant the obs/state/pseudo vectors from hist to forecast
        if cond_type in [:semi, :full]
            T = DSGE.subtract_quarters(date_forecast_start(m), date_prezlb_start(m))

            forecast_output[:histstates] = [x[:, 1:T] for x in histstates]
            forecast_output[:histshocks] = [x[:, 1:T] for x in histshocks]
            forecast_output[:histpseudo] = [x[:, 1:T] for x in histpseudo]
        else
            forecast_output[:histstates] = histstates
            forecast_output[:histshocks] = histshocks
            forecast_output[:histpseudo] = histpseudo
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
            forecast(m, systems, states)

        # For conditional data, transplant the obs/state/pseudo vectors from hist to forecast
        if cond_type in [:semi, :full]
            T = DSGE.subtract_quarters(date_forecast_start(m), date_prezlb_start(m))
            histobs_cond = df_to_matrix(m, df; cond_type = cond_type)[:, index_prezlb_start(m)+T:end]
            
            forecast_output[:forecaststates] = [hcat(x[:, T+1:end], y) for (x, y) in zip(histstates, forecaststates)]
            forecast_output[:forecastshocks] = [hcat(x[:, T+1:end], y) for (x, y) in zip(histshocks, forecastshocks)]
            forecast_output[:forecastpseudo] = [hcat(x[:, T+1:end], y) for (x, y) in zip(histpseudo, forecastpseudo)]
            forecast_output[:forecastobs]    = [hcat(histobs_cond, y) for y in forecastobs]
        else
            forecast_output[:forecaststates] = forecaststates
            forecast_output[:forecastshocks] = forecastshocks
            forecast_output[:forecastpseudo] = forecastpseudo
            forecast_output[:forecastobs]    = forecastobs
        end

        write_forecast_outputs(forecast_vars)
    end


    ### 4. Shock Decompositions

    if !isempty(intersect(output_vars, shockdec_vars))
        shockdecstates, shockdecobs, shockdecpseudo = shock_decompositions(m, systems, histshocks)

        forecast_output[:shockdecstates] = shockdecstates
        forecast_output[:shockdecpseudo] = shockdecpseudo
        forecast_output[:shockdecobs]    = shockdecobs

        write_forecast_outputs(shockdec_vars)
    end


    # Return only saved elements of dict
    filter!((k, v) -> k ∈ output_vars, forecast_output)
    return forecast_output
end


function get_input_file(m, input_type)
    overrides = forecast_input_file_overrides(m)
    if haskey(overrides, input_type)
        override_file = overrides[input_type]
        if ispath(override_file)
            return override_file
        else
            error("Invalid input file override for input_type = $input_type: $override_file")
        end
    end

    if input_type == :mode
        return rawpath(m,"estimate","paramsmode.h5")
    elseif input_type == :mean
        return workpath(m,"estimate","paramsmean.h5")
    elseif input_type == :init
        return ""
    elseif input_type == :full
        return rawpath(m,"estimate","mhsave.h5")
    elseif input_type == :subset
        throw(ArgumentError("Not implemented."))
    else
        throw(ArgumentError("Invalid input_type: $(input_type)"))
    end
end


function get_output_vars(m, output_type)
    if output_type == :states
        vars = [:histstates,
                :histpseudo]
    elseif output_type == :shocks
        vars = [:histshocks]
    elseif output_type == :shocks_nonstandardized
        vars = [:histshocksns]
        throw(ArgumentError("Not implemented."))
    elseif output_type == :forecast
       vars = [:forecaststates,
               :forecastobs,
               :forecastpseudo,
               :forecastshocks]
    elseif output_type == :shockdec
        vars = [:shockdecstates,
                :shockdecpseudo,
                :shockdecobs]
    elseif output_type == :dettrend
        vars = [:dettrendstates,
                :dettrendobs]
        throw(ArgumentError("Not implemented."))
    elseif output_type == :counter
        vars = [:counterstates,
                :counterobs]
        throw(ArgumentError("Not implemented."))
    elseif output_type == :simple
        vars = [:histstates,
                :histpseudo,
                :forecaststates,
                :forecastpseudo,
                :forecastobs,
                :forecastshocks]
    elseif output_type == :all
        vars = []
        throw(ArgumentError("Not implemented."))
    else
        throw(ArgumentError("Invalid input_type: $(output_type)"))
    end
end


function get_output_files(m, input_type, output_vars, cond_type)
    additional_file_strings = ASCIIString[]
    push!(additional_file_strings, "para=" * abbrev_symbol(input_type))
    push!(additional_file_strings, "cond=" * abbrev_symbol(cond_type))

    return [symbol(x) => rawpath(m, "forecast", "$x.jld", additional_file_strings) for x in output_vars]
end
