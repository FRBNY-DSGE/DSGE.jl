"""
```
abstract Scenario
```

Abstract supertype for all alternative scenarios.

### Required Fields

- `key::Symbol`
- `description::String`
- `target_names::Vector{Symbol}`
- `instrument_names::Vector{Symbol}`
- `targets::DataFrame`
- `instruments::DataFrame`
"""
abstract AbstractScenario

type Scenario <: AbstractScenario
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

n_targets(scen::Scenario) = length(scen.target_names)
n_instruments(scen::Scenario) = length(scen.instrument_names)
n_target_horizons(scen::Scenario) = size(scen.targets, 1)

function empty_scenario(key::Symbol, description::String,
                        target_names::Vector{Symbol},
                        instrument_names::Vector{Symbol},
                        vintage::String)
    targets = DataFrame()
    instruments = DataFrame()
    return Scenario(key, description, target_names, instrument_names,
                       targets, instruments, vintage)
end

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

type SwitchingScenario <: AbstractScenario
    key::Symbol
    description::String
    default_scenario_key::Symbol
    prob_enter::Vector{Float64}
    prob_exit::Vector{Float64}
    vintage::String
end

function Base.show(io::IO, scen::SwitchingScenario)
    @printf io "%-12s %s\n" "Original:" scen.key
    @printf io "%-12s %s\n" "Description:" scen.description
    @printf io "%-12s %s\n" "Default:" scen.default_scenario_key
    @printf io "%-12s %s"   "Vintage:" scen.vintage
end

"""
'''
construct_switching_scenario(original, default, prob_enter, prob_exit)
'''
Constructs an instance of `SwitchingScenario` from two scenarios and
two vectors of entry/exit probabilities.
"""
function construct_switching_scenario(original::Scenario, default::Scenario,
                                      prob_enter::Vector{Float64},
                                      prob_exit::Vector{Float64})

    @assert n_target_horizons(original) == length(prob_enter) == length(prob_exit)

    new_scenario = SwitchingScenario(original.key, original.description, default.key,
                                     prob_enter, prob_exit, original.vintage)

    return new_scenario
end
