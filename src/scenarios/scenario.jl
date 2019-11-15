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
abstract type AbstractScenario end

abstract type SingleScenario <: AbstractScenario end

mutable struct Scenario <: SingleScenario
    key::Symbol
    description::String
    target_names::Vector{Symbol}
    instrument_names::Vector{Symbol}
    targets::DataFrame
    instruments::DataFrame
    shock_scaling::Float64
    draw_states::Bool
    altpolicy::AltPolicy
    vintage::String
    n_draws::Int64
end

function Base.show(io::IO, scen::Scenario)
    @printf io "%-14s %s\n" "Key:" scen.key
    @printf io "%-14s %s\n" "Description:" scen.description
    @printf io "%-14s %s\n" "Targets:" scen.target_names
    @printf io "%-14s %s\n" "Instruments:" scen.instrument_names
    if scen.shock_scaling != 1.0
        @printf io "%-14s %s\n" "Shock Scaling:" scen.shock_scaling
    end
    if scen.draw_states
        @printf io "%-14s %s\n" "Draw States:" scen.draw_states
    end
    if scen.altpolicy.key != :historical
        @printf io "%-14s %s\n" "Alt Policy:" scen.altpolicy
    end
    @printf io "%-14s %s\n" "Number of draws:" scen.n_draws
    @printf io "%-14s %s"   "Vintage:" scen.vintage
end

"""
```
Scenario(key, description, target_names, instrument_names, vintage;
    shock_scaling = 1.0, draw_states = false,
    altpolicy = AltPolicy(:historical, eqcond, solve))
```

Scenario constructor. If `instrument_names` is empty, then all model shocks will
be used.
"""
function Scenario(key::Symbol, description::String,
                  target_names::Vector{Symbol},
                  instrument_names::Vector{Symbol},
                  vintage::String;
                  shock_scaling::Float64 = 1.0,
                  draw_states::Bool = false,
                  altpolicy::AltPolicy = AltPolicy(:historical, eqcond, solve),
                  n_draws::Int64 = 0)
    Scenario(key, description, target_names, instrument_names,
             DataFrame(), DataFrame(), shock_scaling, draw_states,
             altpolicy, vintage, n_draws)
end

n_targets(scen::Scenario) = length(scen.target_names)
n_instruments(scen::Scenario) = length(scen.instrument_names)
n_target_horizons(scen::Scenario) = size(scen.targets, 1)
n_scenario_draws(scen::Scenario) = scen.n_draws

function targets_to_data(m::AbstractDSGEModel, scen::Scenario)
    df = DataFrame()

    # Assign dates
    horizons   = n_target_horizons(scen)
    start_date = date_forecast_start(m)
    end_date   = iterate_quarters(start_date, horizons - 1)
    df[!, :date]  = quarter_range(start_date, end_date)

    for var in keys(m.observables)
        df[!, var] = if var in scen.target_names
            scen.targets[!, var]
        else
            fill(missing, horizons)
        end
    end
    return df
end

mutable struct SwitchingScenario <: SingleScenario
    key::Symbol
    description::String
    vintage::String
    original::SingleScenario
    default::SingleScenario
    probs_enter::Vector{Float64}
    probs_exit::Vector{Float64}

    function SwitchingScenario(key, description, vintage, original, default,
                               probs_enter, probs_exit)
        if length(probs_enter) != length(probs_exit)
            error("Lengths of probs_enter and probs_exit must be the same")
        end
        if !all(p -> 0.0 <= p <= 1.0, probs_enter)
            error("Elements of probs_enter must be between 0 and 1")
        end
        if !all(p -> 0.0 <= p <= 1.0, probs_exit)
            error("Elements of probs_exit must be between 0 and 1")
        end
        return new(key, description, vintage, original, default, probs_enter, probs_exit)
    end
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

    SwitchingScenario(key, description, vintage, original, default,
                      probs_enter, probs_exit)
end

"""
```
ScenarioAggregate(key, description, scenarios, sample, proportions, total_draws,
    replace, vintage)

ScenarioAggregate(key, description, scenarios, vintage)

ScenarioAggregate(key, description, scenarios, proportions, total_draws,
    replace, vintage)
```

Composite mutable struct for aggregated `AbstractScenario`s. Scenarios in
`scenarios::Vector{AbstractScenario}` are sampled into the aggregate
with the probability of the corresponding entry in `proportions` (whose values
must sum to 1). The field `replace` indicates whether to sample with
replacement.
"""
mutable struct ScenarioAggregate <: AbstractScenario
    key::Symbol
    description::String
    scenarios::Vector{AbstractScenario}
    sample::Bool # if false, scenario draws are kept in original proportions and `proportions` and `total_draws` are filled in after reading in draws
    proportions::Vector{Float64}
    total_draws::Int
    replace::Bool # whether to sample with replacement
    vintage::String

    function ScenarioAggregate(key, description, scenarios, sample,
                               proportions, total_draws, replace, vintage)
        if sample
            if length(scenarios) != length(proportions)
                error("Lengths of scenarios and proportions must be the same")
            end
            if !all(p -> 0.0 <= p <= 1.0, proportions)
                error("Elements of proportions must be between 0 and 1")
            end
            if !(sum(proportions) ≈ 1.0)
                error("Elements of proportions must sum to 1")
            end
        end
        return new(key, description, scenarios, sample,
                   proportions, total_draws, replace, vintage)
    end
end

# `sample = true` constructor
function ScenarioAggregate(key::Symbol, description::String, scenarios::Vector{AbstractScenario},
                           proportions::Vector{Float64}, total_draws::Int, replace::Bool,
                           vintage::String)
    sample = true
    return ScenarioAggregate(key, description, scenarios, sample, proportions, total_draws,
                             replace, vintage)
end

# `sample = false` constructor
function ScenarioAggregate(key::Symbol, description::String, scenarios::Vector{AbstractScenario},
                           vintage::String)
    sample = false
    proportions = Float64[]
    total_draws = -1
    replace = false
    return ScenarioAggregate(key, description, scenarios, sample, proportions, total_draws,
                             replace, vintage)
end

function Base.show(io::IO, agg::ScenarioAggregate)
    @printf io "%-24s %s\n" "Key:" agg.key
    @printf io "%-24s %s\n" "Description:" agg.description
    @printf io "%-24s %s\n" "Scenarios:" map(scen -> scen.key, agg.scenarios)
    @printf io "%-24s %s\n" "Sample:" agg.sample
    @printf io "%-24s %s\n" "Proportions:" agg.proportions
    @printf io "%-24s %s\n" "Total Draws:" agg.total_draws
    @printf io "%-24s %s\n" "Sample with Replacement:" agg.replace
    @printf io "%-24s %s"   "Vintage:" agg.vintage
end
