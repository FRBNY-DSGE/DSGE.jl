"""
```
compute_scenario_system(m, scen::Scenario; apply_altpolicy = false)
```

Given the current model parameters, compute the state-space system corresponding
to model `m` and alternative scenario `scen`. This function differs from
`compute_system` in that the `CCC`, `DD`, and `DD_pseudo` vectors are set to
zero (since we forecast in deviations from baseline) and shocks that are not in
`scen.instrument_names` are zeroed out in the `QQ` matrix.
"""
function compute_scenario_system(m::AbstractDSGEModel, scen::Scenario;
                                 tvis::Bool = false, apply_altpolicy::Bool = false)

    system = compute_system(m, tvis = tvis)

    # Set constant system matrices to 0
    system = zero_system_constants(system)

    # Zero out non-instrument shocks
    system.measurement.QQ = copy(system[:QQ])
    for shock in keys(m.exogenous_shocks)
        if !(shock in scen.instrument_names)
            shock_index = m.exogenous_shocks[shock]
            system[:QQ][shock_index, :] .= 0
            system[:QQ][:, shock_index] .= 0
        end
    end

    # Error out if there is nonzero measurement error
    if any(x -> x != 0, system[:EE])
        error("Can't simulate scenarios under nonzero measurement error")
    end

    return system
end

"""
```
filter_shocks!(m, scen, system::Scenario)
```

Given a scenario draw `scen`, back out the shocks necessary to hit
`scen.targets` and put them into `scen.instruments`. This function returns
`forecastshocks`, an `nshocks` x `horizon` matrix of filtered and smoothed
shocks.

This function checks `forecast_uncertainty_override(m)` for whether to smooth
shocks using the simulation smoother.
"""
function filter_shocks!(m::AbstractDSGEModel, scen::Scenario, system::System, in_sample::Bool = false)
    # Check applicability of this methodology
    @assert n_instruments(scen) >= n_targets(scen) "Number of instruments must be at least number of targets"

    # Construct data
    df = targets_to_data(m, scen)

    # Set initial state and state covariance to 0
    s_0 = zeros(n_states_augmented(m))
    P_0 = zeros(n_states_augmented(m), n_states_augmented(m))

    # Filter and smooth *deviations from baseline*
    _, forecastshocks, _ = smooth(m, df, system, s_0, P_0, draw_states = scen.draw_states,
                                  include_presample = true, in_sample = in_sample)

    # Assign shocks to instruments DataFrame
    for shock in keys(m.exogenous_shocks)
        shock_index = m.exogenous_shocks[shock]
        if shock in scen.instrument_names
            scen.instruments[!, shock] = forecastshocks[shock_index, :]
        else
            @assert all(x -> x == 0, forecastshocks[shock_index, :])
        end
    end

    return forecastshocks
end

"""
```
forecast_scenario_draw(m, scen::Scenario, system, draw_index)
```

Filter shocks and use them to forecast the `draw_index`th draw of `scen`.
"""
function forecast_scenario_draw(m::AbstractDSGEModel, scen::Scenario, system::System,
                                draw_index::Int)
    # Re-initialize model indices in case extra states or equations were added
    # for an alternative policy
    init_model_indices!(m)

    # Load targets
    load_scenario_targets!(m, scen, draw_index)

    # Filter shocks
    forecastshocks = filter_shocks!(m, scen, system)

    # Scale shocks if desired
    forecastshocks = scen.shock_scaling * forecastshocks

    # Re-solve model with alternative policy rule, if applicable
    if alternative_policy(m).key != :historical
        system = compute_scenario_system(m, scen, apply_altpolicy = true)
    end

    # Forecast
    s_T = zeros(n_states_augmented(m))
    forecaststates, forecastobs, forecastpseudo, _ =
        forecast(m, system, s_T, shocks = forecastshocks)

    # Check forecasted output matches targets *if not forecasting under
    # alternative policy or using simulation smoother*
    if alternative_policy(m).key == :historical && !scen.draw_states
        for var in scen.target_names
            var_index = m.observables[var]
            horizon = min(forecast_horizons(m), n_target_horizons(scen))
            @assert forecastobs[var_index, 1:horizon] ≈ scen.shock_scaling * scen.targets[1:horizon, var]
        end
    end

    # Return output dictionary
    forecast_output = Dict{Symbol, Array{Float64}}()
    forecast_output[:forecaststates] = forecaststates
    forecast_output[:forecastobs]    = forecastobs
    forecast_output[:forecastpseudo] = forecastpseudo
    forecast_output[:forecastshocks] = forecastshocks

    return forecast_output
end

"""
```
write_scenario_forecasts(m, scenario_output_files, forecast_output;
    verbose = :low)
```

Write scenario outputs in `forecast_output` to `values(scenario_output_files)`.
"""
function write_scenario_forecasts(m::AbstractDSGEModel,
                                  scenario_output_files::Dict{Symbol, String},
                                  forecast_output::Dict{Symbol, Array{Float64}};
                                  verbose::Symbol = :low)
    for (i, var) in enumerate([:forecastobs, :forecastpseudo])
        filepath = scenario_output_files[var]
        JLD2.jldopen(filepath, "w") do file
            write_forecast_metadata(m, file, var)
            write(file, "arr", Array{Float64}(forecast_output[var]))
            if :proportion_switched in keys(scenario_output_files)
                write(file, "proportion_switched", forecast_output[:proportion_switched][i])
            end
        end

        println(verbose, :high, " * Wrote " * basename(filepath))
    end
end

"""
```
forecast_scenario(m, scen::Scenario; verbose = :low)
```

Simulate all draws of `scen` using the modal parameters of the model `m`. This
function returns a `Dict{Symbol, Array{Float64}`.
"""
function forecast_scenario(m::AbstractDSGEModel, scen::Scenario;
                           verbose::Symbol = :low)
    # Print
    info_print(verbose, :low, "Forecasting scenario = " * string(scen.key) * "...")
    println(verbose, :low, "Start time: " * string(now()))
    println(verbose, :low, "Forecast outputs will be saved in " * rawpath(m, "scenarios"))

    start_time = time_ns()
    # Update model alt policy setting
    m <= Setting(:alternative_policy, scen.altpolicy, false, "apol",
                 "Alternative policy")

    # If no instrument names provided, use all shocks
    if isempty(scen.instrument_names)
        scen.instrument_names = collect(keys(m.exogenous_shocks))
    end

    # Load modal parameters and compute system
    params = load_draws(m, :mode; verbose = verbose)
    update!(m, params)
    system = compute_scenario_system(m, scen)

    # Get to work!
    ndraws = scen.n_draws == 0 ? count_scenario_draws!(m, scen) : scen.n_draws
    mapfcn = use_parallel_workers(m) ? pmap : map
    forecast_outputs = mapfcn(draw_ind -> forecast_scenario_draw(m, scen, system, draw_ind),
                              1:ndraws)

    # Assemble outputs and write to file
    forecast_outputs = convert(Vector{Dict{Symbol, Array{Float64}}}, forecast_outputs)
    forecast_output = assemble_block_outputs(forecast_outputs)
    output_files = get_scenario_output_files(m, scen, [:forecastobs, :forecastpseudo])
    write_scenario_forecasts(m, output_files, forecast_output, verbose = verbose)
    # Print
    forecast_time = (time_ns() - start_time)/1e9
    forecast_time_min = forecast_time/60
    println(verbose, :low, "\nTime elapsed: " * string(forecast_time_min) * " minutes")
    println(verbose, :low, "Forecast complete: " * string(now()))

    return forecast_output
end
