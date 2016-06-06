"""
```
forecast_all(m::AbstractModel, df::DataFrame; cond_types::Vector{Symbol}
input_types::Vector{Symbol} output_types::Vector{Symbol}
```

Compute forecasts for all specified combinations of conditional data, input types, and
output types.

# Arguments

- `m`: model object
- `df`: DataFrame of data for observables
- `cond_types`: conditional data type, any combination of
    - `:none`: no conditional data
    - `:semi`: use "semiconditional data" - average of quarter-to-date observations for high frequency series
    - `:full`: use "conditional data" - semiconditional plus nowcasts for desired
      observables
- `input_types`: which set of parameters to use, any combination of
    - `:mode`: forecast using the modal parameters only
    - `:mean`: forecast using the mean parameters only
    - `:init`: forecast using the initial parameter values only
    - `:full`: forecast using all parameters (full distribution)
    - `:subset`: forecast using a well-defined user-specified subset of draws
- `output_types`: forecast routine outputs to compute, any combination of
    - `:states`: smoothed states (history) for all specified conditional data types
    - `:shocks`: smoothed shocks (history, standardized) for all specified conditional data types
    - `:shocks_nonstandardized`: smoothed shocks (history, non-standardized) for all
        specified conditional data types
    - `:forecast`: forecast of states and observables for all specified conditional data types
    - `:shockdec`: shock decompositions (history) of states and observables for all
        specified conditional data types
    - `:dettrend`: deterministic trend (history) of states and observables for all specified
        conditional data types
    - `:counter`: counterfactuals (history) of states and observables for all specified
        conditional data types
    - `:simple`: smoothed states, forecast of states, forecast of observables for
        *unconditional* data only
    - `:simple_cond`: smoothed states, forecast of states, forecast of observables for all
        specified conditional data types
    - `:all`: smoothed states (history), smoothed shocks (history, standardized), smoothed
      shocks (history, non-standardized), shock decompositions (history), deterministic
      trend (history), counterfactuals (history), forecast, forecast shocks drawn, shock
      decompositions (forecast), deterministic trend (forecast), counterfactuals (forecast)
   Note that some similar outputs may or may not fall under the "forecast_all" framework,
   including
    - `:mats`: recompute system matrices (TTT, RRR, CCC) given parameters only
    - `:zend`: recompute final state vector (s_{T}) given parameters only
    - `:irfs`: impulse response functions

Outputs
-------

- todo
"""
function forecast_all(m::AbstractModel, df::DataFrame;
                      cond_types::Vector{Symbol}   = Vector{Symbol}(),
                      input_types::Vector{Symbol}  = Vector{Symbol}(),
                      output_types::Vector{Symbol} = Vector{Symbol}())

    for input_type in input_types
        for output_type in output_types
            for cond_type in cond_types
                forecast_one(m, df; cond_type=cond_type, input_type=input_type, output_type=output_type)
            end
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
prepare_states(m::AbstractModel, input_type::Symbol, systems::Vector{System{Float64}},
params::Matrix{Float64}, df::DataFrame, zend::Matrix{Float64})
```
"""
function prepare_states(m::AbstractModel, input_type::Symbol,
    systems::Vector{System{Float64}}, params::Matrix{Float64}, df::DataFrame,
    zend::Matrix{Float64})

    # Setup and preallocate
    n_sim_forecast = size(systems,1)
    n_sim = size(params,1)
    jstep = convert(Int, n_sim/n_sim_forecast)
    states = Vector{Vector{Float64}}(n_sim_forecast)

    # If we just have one draw of parameters in mode, mean, or init case, then we don't have the
    # pre-computed system matrices. We now recompute them here by running the Kalman filter.
    if input_type in [:mean, :mode, :init]
        update!(m, vec(params))
        filt, _, _ = filter(m, df, systems; Ny0 = n_presample_periods(m))
        # filt is a vector of nperiods x nstates matrices of filtered states
        states[1] = vec(filt[1][end,:])

    # If we have many draws, then we must package them into a vector of System objects.
    elseif input_type in [:full]
        for i in 1:n_sim_forecast
            j = i * jstep
            states[i] = vec(zend[j,:])
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
        systems[1] = compute_system(m; use_expected_rate_data = true)
    elseif input_type in [:full]
        empty = isempty(CCC)
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
prepare_forecast_inputs(m::AbstractModel, df::DataFrame; input_type::Symbol  = :mode,
output_type::Symbol = :simple, cond_type::Symbol  = :none)
```

Load draws for this input type, prepare a System object for each draw, and prepare initial
state vectors.
"""
function prepare_forecast_inputs(m::AbstractModel, df::DataFrame;
                      input_type::Symbol  = :mode,
                      output_type::Symbol = :simple,
                      cond_type::Symbol  = :none)
    # Some variables
    n_states = n_states_augmented(m)

    # Set up infiles
    println("Loading draws")
    params, TTT, RRR, CCC, zend = load_draws(m, input_type)

    n_sim = size(params,1)
    jstep = get_jstep(m, n_sim)
    n_sim_forecast = convert(Int, n_sim/jstep)

    # Populate systems vector
    println("Preparing systems")
    systems = prepare_systems(m, input_type, params, TTT, RRR, CCC)

    # Populate states vector
    println("Preparing states")
    states = prepare_states(m, input_type, systems, params, df, zend)

    return systems, states
end

function forecast_one(m::AbstractModel, df::DataFrame;
                      input_type::Symbol  = :mode,
                      output_type::Symbol = :simple,
                      cond_type::Symbol  = :none)

    systems, states = prepare_forecast_inputs(m, df;
        input_type=input_type, output_type=output_type, cond_type=cond_type)

    # Prepare conditional data matrix. All missing columns will be set to NaN.
    if cond_type in [:semi, :full]
        cond_data = load_cond_data(m, cond_type)
        cond_df = [df; cond_data]
    end

    # Example: call forecast, unconditional data, states+observables
    forecast_output = Dict{Symbol, Vector{Array{Float64}}}()

    if output_type in [:forecast, :simple, :simple_cond]
        println("Calling forecast")
        forecaststates, forecastobs, forecastpseudo = 
            forecast(m, systems, states)

        forecast_output[:forecastobs] = forecastobs
        forecast_output[:forecaststates] = forecaststates
        forecast_output[:forecastpseudo] = forecastpseudo
        # forecast_output[:forecastshocks] = forecastshocks

    end

    # Set up outfiles
    output_files = get_output_files(m, input_type, output_type, cond_type)

    # Write output files
    println("Writing output files")
    for (var,file) in output_files
        jldopen(file, "w") do f
            write(f, string(var), forecast_output[var])
        end
    end

    return forecast_output

end

    
function get_input_file(m, input_type)
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

function get_output_files(m, input_type, output_type, cond_type)

    # Add additional file strings here
    additional_file_strings = ASCIIString[]
    push!(additional_file_strings, "para=" * abbrev_symbol(input_type))
    push!(additional_file_strings, "cond=" * abbrev_symbol(cond_type))

    # vars prefix
    if output_type == :states
        vars = ["histstates"]
        throw(ArgumentError("Not implemented."))
    elseif output_type == :shocks
        vars = ["histshocks"]
        throw(ArgumentError("Not implemented."))
    elseif output_type == :shocks_nonstandardized
        vars = ["histshocksns"]
        throw(ArgumentError("Not implemented."))
    elseif output_type == :forecast
       vars = ["forecaststates",
               "forecastobs",
               "forecastpseudo"]
#                   "forecastshocks"]
    elseif output_type == :shockdec
        vars = ["shockdecstates",
                   "shockdecobs"]
        throw(ArgumentError("Not implemented."))
    elseif output_type == :dettrend
        vars = ["dettrendstates",
                   "dettrendobs"]
        throw(ArgumentError("Not implemented."))
    elseif output_type == :counter
        vars = ["counterstates",
                   "counterobs"]
        throw(ArgumentError("Not implemented."))
    elseif output_type in [:simple, :simple_cond]
        vars = ["histstates",
                   "forecaststates",
                   "forecastobs",
                   "forecastshocks"]
        throw(ArgumentError("Not implemented."))
    elseif output_type == :all
        vars = []
        throw(ArgumentError("Not implemented."))
    else
        throw(ArgumentError("Invalid input_type: $(output_type)"))
    end

    return [symbol(x) => rawpath(m, "forecast", x*".jld", additional_file_strings) for x in vars]
end
