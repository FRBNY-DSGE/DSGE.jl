abstract Scenario

function Base.show(io::IO, scen::Scenario)
    @printf io "Scenario: %s\n" scen.description
    @printf io "Targets: %s\n" scen.target_names
    @printf io "Instruments: %s\n" scen.instrument_names
end

n_targets(scen::Scenario) = length(scen.target_names)
n_instruments(scen::Scenario) = length(scen.instrument_names)
n_target_horizons(scen::Scenario) = size(scen.targets, 1)

function scenario_targets_file(m::AbstractModel)
    key = get_setting(m, :scenario_key)
    vint = get_setting(m, :scenario_vintage)
    basename = string(key) * "_" * string(vint) * ".jld"
    return inpath(m, "scenarios", basename)
end

function load_scenario_targets!(scen::Scenario, path::String, draw_index::Int)
    raw_targets = h5read(path, "arr", (:, :, draw_index))
    target_inds = load(path, "target_indices")

    @assert collect(keys(target_inds)) == scen.target_names "Target indices in $path do not match target names in $(scen.key)"

    for (target_name, target_index) in target_inds
        scen.targets[target_name] = raw_targets[target_index, :]
    end

    return scen.targets
end

function targets_to_data(m::AbstractModel, scen::Scenario)
    df = DataFrame()
    horizons = n_target_horizons(scen)
    df[:date] = DSGE.quarter_range(date_forecast_start(m), DSGE.iterate_quarters(date_forecast_start(m), horizons - 1))
    for var in keys(m.observables)
        df[var] = if var in scen.target_names
            scen.targets[var]
        else
            fill(NaN, horizons)
        end
    end
    return df
end

type MyScenario <: Scenario
    key::Symbol
    description::String
    target_names::Vector{Symbol}
    instrument_names::Vector{Symbol}
    targets::DataFrame
    instruments::DataFrame
end

function MyScenario()
    key = :myscenario
    description = "My Scenario"
    target_names = [:obs_gdp, :obs_corepce]
    instrument_names = [:g_sh, :b_sh, :μ_sh, :z_sh]
    targets = DataFrame()
    instruments = DataFrame()
    return MyScenario(key, description, target_names, instrument_names, targets, instruments)
end