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

function forecast_one(m::AbstractModel, df::DataFrame;
                      input_type::Symbol  = :mode,
                      output_type::Symbol = :simple,
                      cond_type::Symbol  = :none)
    # Some variables
    n_states = n_states_augmented(m)
    jstep = get_setting(m, :forecast_jstep)

    # Set up infiles
    input_file_name = get_input_file(m, input_type)

    # Read infiles and set n_sim based on input_type type
    if input_type in [:mean, :mode]
        h5open(input_file_name, "r") do f
            params = map(Float64, read(f, "params"))
        end
        TTT = RRR = CCC = zend = []
        n_sim = 1
        jstep = 1
    elseif input_type == :full
        params, TTT, RRR, zend = h5open(input_file_name, "r") do f
            params = map(Float64, read(f, "mhparams"))
            TTT    = map(Float64, read(f, "mhTTT"))
            RRR    = map(Float64, read(f, "mhRRR"))
            #CCC   = map(Float64, read(f, "mhCCC"))
            zend   = map(Float64, read(f, "mhzend"))
            params, TTT, RRR, zend
        end
        n_sim = size(params,1)
    end

    n_sim_forecast = convert(Int, n_sim/jstep)

    # Populate systems vector
    systems = Vector{System{Float64}}(n_sim_forecast)
    initial_state_draws = Vector{Vector{Float64}}(n_sim_forecast)

    # If we just have one draw of parameters in mode or mean case, then we don't have the
    # pre-computed system matrices. We now recompute them here by running the Kalman filter.
    if input_type in [:mean, :mode]
        update!(m, params)
        sys = compute_system(m; use_expected_rate_data = true)
        kal = filter(m, df, sys; Ny0 = n_presample_periods(m))
        zend = kal[:zend]

        # Prepare system
        systems[1] = sys
        initial_state_draws[1] = vec(zend)

    # If we have many draws, then we must package them into a vector of System objects.
    elseif input_type in [:full]
        for i in 1:n_sim_forecast
            j = i * jstep
            # Prepare transition eq
            TTT_j  = squeeze(TTT[j,:,:],1)
            RRR_j  = squeeze(RRR[j,:,:],1)
            #CCC_j = squeeze(CCC[j,:,:],1)
            trans_j = Transition(TTT_j, RRR_j)
            CCC_j = trans_j[:CCC]

            # Prepare measurement eq
            params_j = vec(params[j,:])
            update!(m, params_j)
            meas_j   = measurement(m, TTT_j, RRR_j, CCC_j; shocks = true)

            # Prepare system
            sys_j = System(trans_j, meas_j)
            systems[i] = sys_j
            initial_state_draws[i] = vec(zend[j,:])
        end
    end

    # Prepare conditional data matrix. All missing columns will be set to NaN.
    if cond_type in [:semi, :full]
        cond_data = load_cond_data(m, cond_type)
        df = [df; cond_data]
    end

    # Example: call forecast, unconditional data, states+observables
    forecast_output = Dict{Symbol, Any}()

    if output_type in [:forecast, :simple, :simple_cond]
        forecastobs, forecaststates, forecastshocks = forecast(m, systems, initial_state_draws)
        forecast_output[:forecastobs] = forecastobs
        forecast_output[:forecaststates] = forecaststates
        forecast_output[:forecastshocks] = forecastshocks
    end

    # Set up outfiles
    output_files = get_output_files(m, input_type, output_type, cond_type)

    # Write output files
    for (var,file) in output_files
        h5open(file, "w") do f
            write(f, string(var), forecast_output[var])
        end
    end

end

function get_input_file(m, input_type)
    if input_type == :mode
        return rawpath(m,"estimate","paramsmode.h5")
    elseif input_type == :mean
        return workpath(m,"estimate","paramsmean.h5")
    elseif input_type == :full
        return rawpath(m,"estimate","mhsave.h5")
    elseif input_type == :subset
        #TODO
        return ""
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
                   "forecastshocks"]
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
    end

    return [symbol(x) => rawpath(m, "forecast", x*".h5", additional_file_strings) for x in vars]
end
