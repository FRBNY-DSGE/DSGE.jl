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
    systems = Vector{System}(n_sim_forecast)

    # If we just have one draw of parameters in mode or mean case, then we don't have the
    # pre-computed system matrices. We now recompute them here by running the Kalman filter.
    if input_type in [:mean, :mode]
        update!(m, params)
        sys = compute_system(m; use_expected_rate_date = false)
        kal = filter(m, df, sys; Ny0 = n_presample_periods(m))
        zend = kal[:zend]
        initial_state_draws = vcat(zend)

    # If we have many draws, then we must package them into a vector of System objects.
    elseif input_type in [:full]
        for i in 1:n_sim_forecast
            j = i * jstep
            # Prepare transition eq
            TTT_j  = squeeze(TTT[j,:,:],1)
            RRR_j  = squeeze(RRR[j,:,:],1)
            #CCC_j = squeeze(CCC[j,:,:],1)
            CCC_j  = zeros(eltype(TTT_j), size(TTT_j, 1), 1)
            trans_j = Transition(TTT_j, RRR_j, CCC_j)

            # Prepare measurement eq
            params_j = vec(params[j,:])
            update!(m, params_j)
            meas_j   = measurement(m, TTT_j, RRR_j, CCC_j; shocks = false)

            # Prepare system
            sys_j = System(trans_j, meas_j)
            systems[i] = sys_j
        end

        initial_state_draws = zend[jstep:jstep:n_sim,:]
    end

    # Check initial_state_draws matrix, which is n_simulations x n_states
    @assert size(initial_state_draws) == (n_sim/jstep, n_states)

    # Prepare conditional data matrix. All missing columns will be set to NaN.
    if cond_type in [:semi, :full]
        cond_data = load_cond_data(m, cond_type)
        df = [df; cond_data]
    end

    # Example: call forecast, unconditional data, states+observables
    function forecast(args...)
        return 1,2,3
    end
    forecastobs, forecaststates, forecastshocks = forecast(m, systems, initial_state_draws)

    # Set up outfiles
    output_file_names = get_output_files(m, input_type, output_type, cond_type)

    # In this demo, output_file_names is a three element vector
    h5open(output_file_names[1], "w") do f
        f["forecastobs"] = forecastobs
    end
    h5open(output_file_names[2], "w") do f
        f["forecaststates"] = forecaststates
    end
    h5open(output_file_names[3], "w") do f
        f["forecastshocks"] = forecastshocks
    end

    # # Write outfiles
    # for output_file in output_file_names
    #     write(output_file)
    # end

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

function get_output_files(m, input_type, output_type, cond)

    # Add additional file strings here
    additional_file_strings = []
    push!(additional_file_strings, "para=" * abbrev_symbol(input_type))
    push!(additional_file_strings, "cond=" * abbrev_symbol(cond))

    # Results prefix
    if output_type == :states
        results = ["histstates"]
        throw(ArgumentError("Not implemented."))
    elseif output_type == :shocks
        results = ["histshocks"]
        throw(ArgumentError("Not implemented."))
    elseif output_type == :shocks_nonstandardized
        results = ["histshocksns"]
        throw(ArgumentError("Not implemented."))
    elseif output_type == :forecast
        results = ["forecaststates",
                   "forecastobs",
                   "forecastshocks"]
    elseif output_type == :shockdec
        results = ["shockdecstates",
                   "shockdecobs"]
        throw(ArgumentError("Not implemented."))
    elseif output_type == :dettrend
        results = ["dettrendstates",
                   "dettrendobs"]
        throw(ArgumentError("Not implemented."))
    elseif output_type == :counter
        results = ["counterstates",
                   "counterobs"]
        throw(ArgumentError("Not implemented."))
    elseif output_type in [:simple, :simple_cond]
        results = ["histstates",
                   "forecaststates",
                   "forecastobs",
                   "forecastshocks"]
        throw(ArgumentError("Not implemented."))
    elseif output_type == :all
        results = []
        throw(ArgumentError("Not implemented."))
    end

    return [rawpath(m, "forecast", x*".h5", additional_file_strings) for x in results]
end
