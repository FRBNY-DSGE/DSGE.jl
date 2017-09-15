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
abstract Scenario

function Base.show(io::IO, scen::Scenario)
    @printf io "%-12s %s\n" "Key:" scen.key
    @printf io "%-12s %s\n" "Description:" scen.description
    @printf io "%-12s %s\n" "Targets:" scen.target_names
    @printf io "%-12s %s" "Instruments:" scen.instrument_names
end

n_targets(scen::Scenario) = length(scen.target_names)
n_instruments(scen::Scenario) = length(scen.instrument_names)
n_target_horizons(scen::Scenario) = size(scen.targets, 1)

function empty_scenario(constructor::DataType, key::Symbol, description::String,
                        target_names::Vector{Symbol}, instrument_names::Vector{Symbol})
    targets = DataFrame()
    instruments = DataFrame()
    return constructor(key, description, target_names, instrument_names,
                       targets, instruments)
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

abstract SwitchingScenario <: Scenario

function Base.show(io::IO, scen::SwitchingScenario)
    @printf io "%-12s %s\n" "Original:" scen.key
    @printf io "%-12s %s\n" "Description:" scen.description
    @printf io "%-12s %s\n" "Description:" scen.default_scenario
end


function define_switching_scenario(original::Scenario)
    new_name = Symbol(string(original.key) * "Switching")

    return @eval begin
        type $(new_name) <: SwitchingScenario
            key::Symbol
            description::String
            default_scenario::Symbol
            probs_enter::Vector{Float64}
            probs_leave::Vector{Float64}
        end
    end
 end

function construct_switching_scenario(original::Scenario, default::Scenario,
                                      prob_enter::Vector{Float64},
                                      prob_leave::Vector{Float64})

    new_name = Symbol(string(original.key) * "Switching")
    if !isdefined(new_name)
        define_switching_scenario(original)
    end

    @assert length(prob_enter) == length(prob_leave) == n_target_horizons(original)

    @eval constructor = $(new_name)
    new_scenario = constructor(original.key, original.description,
                               default.key, prob_enter, prob_leave)
    return new_scenario
end
