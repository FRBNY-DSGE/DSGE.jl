function compute_scenario_system(m::AbstractModel, scen::Scenario)
    system = compute_system(m)

    # Set D = 0
    system.measurement.DD = zeros(size(system[:DD]))

    # Zero out non-instrument shocks
    system.measurement.QQ = copy(system[:QQ])
    for shock in keys(m.exogenous_shocks)
        if !(shock in scen.instrument_names)
            shock_index = m.exogenous_shocks[shock]
            system[:QQ][shock_index, :] = 0
            system[:QQ][:, shock_index] = 0
        end
    end

    return system
end

function filter_shocks!(m::AbstractModel, scen::Scenario, system::System)
    # Check applicability of this methodology
    @assert n_instruments(scen) >= n_targets(scen) "Number of instruments must be at least number of targets"

    # Construct data
    df = targets_to_data(m, scen)

    # Set initial state and state covariance to 0
    s_0 = zeros(n_states_augmented(m))
    P_0 = zeros(n_states_augmented(m), n_states_augmented(m))

    # Filter and smooth *deviations from baseline*
    kal = DSGE.filter(m, df, system, s_0, P_0)
    _, forecastshocks, _ = smooth(m, df, system, kal, draw_states = false,
                                  include_presample = true)

    # Assign shocks to instruments DataFrame
    for shock in keys(m.exogenous_shocks)
        shock_index = m.exogenous_shocks[shock]
        if shock in scen.instrument_names
            scen.instruments[shock] = forecastshocks[shock_index, :]
        else
            @assert all(x -> x == 0, forecastshocks[shock_index, :])
        end
    end

    return forecastshocks
end

function forecast_scenario_draw(m::AbstractModel, scenario_key::Symbol, scenario_vint::String,
                                draw_index::Int)
    # Initialize scenario
    constructor = eval(scenario_key)
    scen = constructor()
    load_scenario_targets!(m, scen, scenario_vint, draw_index)

    # Filter shocks
    # TODO: add solving for shocks
    system = compute_scenario_system(m, scen)
    forecastshocks = filter_shocks!(m, scen, system)

    # Forecast
    s_T = zeros(n_states_augmented(m))
    forecaststates, forecastobs, forecastpseudo, _ =
        forecast(m, system, s_T, shocks = forecastshocks)

    # Return output dictionary
    forecast_output = Dict{Symbol, Array{Float64}}()
    forecast_output[:forecaststates] = forecaststates
    forecast_output[:forecastobs]    = forecastobs
    forecast_output[:forecastpseudo] = forecastpseudo
    forecast_output[:forecastshocks] = forecastshocks

    return forecast_output
end

function write_scenario_forecasts(m::AbstractModel,
                                  scenario_output_files::Dict{Symbol, String},
                                  forecast_output::Dict{Symbol, Array{Float64}};
                                  verbose::Symbol = :low)
    for var in [:forecastobs, :forecastpseudo]
        filepath = scenario_output_files[var]
        jldopen(filepath, "w") do file
            DSGE.write_forecast_metadata(m, file, var)
            write(file, "arr", forecast_output[var])
        end

        if DSGE.VERBOSITY[verbose] >= DSGE.VERBOSITY[:high]
            println(" * Wrote $(basename(filepath))")
        end
    end
end

function forecast_scenario(m::AbstractModel, scenario_key::Symbol,
                           scenario_vint::String; verbose::Symbol = :low)
    # Print
    if DSGE.VERBOSITY[verbose] >= DSGE.VERBOSITY[:low]
        info("Forecasting scenario = $scenario_key...")
        println("Start time: " * string(now()))
        println("Forecast outputs will be saved in " * rawpath(m, "scenarios"))
    end

    # Load modal parameters
    params = load_draws(m, :mode; verbose = verbose)
    DSGE.update!(m, params)

    # Get to work!
    ndraws = n_scenario_draws(m, scenario_key, scenario_vint)
    mapfcn = use_parallel_workers(m) ? pmap : map
    forecast_outputs = mapfcn(draw_ind -> forecast_scenario_draw(m, scenario_key, scenario_vint, draw_ind),
                              1:ndraws)

    # Assemble outputs and write to file
    forecast_outputs = convert(Vector{Dict{Symbol, Array{Float64}}}, forecast_outputs)
    forecast_output = DSGE.assemble_block_outputs(forecast_outputs)
    output_files = get_scenario_output_files(m, scenario_key, scenario_vint, [:forecastobs, :forecastpseudo])
    write_scenario_forecasts(m, output_files, forecast_output, verbose = verbose)

    return forecast_output
end