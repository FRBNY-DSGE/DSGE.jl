"""
```
abstract AbstractScenario
```

Abstract supertype for all alternative scenarios.

### Required Fields

- `key::Symbol`
- `description::String`
- `vintage::String`
"""
abstract AbstractScenario

abstract SingleScenario <: AbstractScenario

type Scenario <: SingleScenario
    key::Symbol
    description::String
    target_names::Vector{Symbol}
    instrument_names::Vector{Symbol}
    targets::DataFrame
    instruments::DataFrame
    vintage::String
end

function Base.show(io::IO, scen::Scenario)
    @printf io "%-12s %s\n" "Key:" scen.key
    @printf io "%-12s %s\n" "Description:" scen.description
    @printf io "%-12s %s\n" "Targets:" scen.target_names
    @printf io "%-12s %s\n" "Instruments:" scen.instrument_names
    @printf io "%-12s %s"   "Vintage:" scen.vintage
end

function Scenario(key::Symbol, description::String,
                  target_names::Vector{Symbol},
                  instrument_names::Vector{Symbol},
                  vintage::String)
    targets = DataFrame()
    instruments = DataFrame()
    return Scenario(key, description, target_names, instrument_names,
                       targets, instruments, vintage)
end

n_targets(scen::Scenario) = length(scen.target_names)
n_instruments(scen::Scenario) = length(scen.instrument_names)
n_target_horizons(scen::Scenario) = size(scen.targets, 1)

function targets_to_data(m::AbstractModel, scen::Scenario)
    df = DataFrame()

    # Assign dates
    horizons = n_target_horizons(scen)
    start_date = date_forecast_start(m)
    end_date = DSGE.iterate_quarters(date_forecast_start(m), horizons - 1)
    df[:date] = DSGE.quarter_range(start_date, end_date)

    for var in keys(m.observables)
        df[var] = if var in scen.target_names
            scen.targets[var]
        else
            fill(NaN, horizons)
        end
    end
    return df
end

type SwitchingScenario <: SingleScenario
    key::Symbol
    description::String
    vintage::String
    original::SingleScenario
    default::SingleScenario
    probs_enter::Vector{Float64}
    probs_exit::Vector{Float64}
end

function Base.show(io::IO, scen::SwitchingScenario)
    @printf io "%-12s %s\n" "Key:" scen.key
    @printf io "%-12s %s\n" "Description:" scen.description
    @printf io "%-12s %s\n" "Original:" scen.original.key
    @printf io "%-12s %s\n" "Default:" scen.default.key
    @printf io "%-12s %s"   "Vintage:" scen.vintage
end

"""
```
SwitchingScenario(key, original, default, probs_enter, probs_exit;
                  description = original.description,
                  vintage = original.vintage)
```

Constructs an instance of `SwitchingScenario` from two scenarios and
two vectors of entry/exit probabilities.
"""
function SwitchingScenario(key::Symbol, original::Scenario, default::Scenario,
                           probs_enter::Vector{Float64},
                           probs_exit::Vector{Float64};
                           description::String = original.description,
                           vintage::String = original.vintage)

    @assert n_target_horizons(original) <= length(probs_enter)
    @assert n_target_horizons(original) <= length(probs_exit)

    new_scenario = SwitchingScenario(key, description, vintage,
                                     original, default,
                                     probs_enter, probs_exit)

    return new_scenario
end

"""
```
type ScenarioAggregate
```

Composite type for aggregated alternative scenarios.
"""
type ScenarioAggregate <: AbstractScenario
    key::Symbol
    description::String
    scenario_groups::Vector{Vector{SingleScenario}}
    proportions::Vector{Float64}
    total_draws::Int
    replace::Bool # whether to sample with replacement
    vintage::String
end
