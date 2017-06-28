"""
```
AltPolicy(rule, [forecast_init = identity])
```

Type defining an alternative policy rule.

### Fields

- `rule::Function`: A function that solves the model and replaces the
  baseline policy rule with the desired alternative rule. It is also
  responsible for augmenting the state-space system with the states in
  `m.endogenous_states_augmented`. This function must accept 1
  argument: an instance of a subtype of `AbstractModel`. It must return

- `forecast_init::Function`: A function that initializes forecasts
  under the alternative policy rule. Specifically, it accepts a model,
  an `nshocks` x `n_forecast_periods` matrix of shocks to be applied
  in the forecast, and a vector of initial states for the forecast. It
  must return a new matrix of shocks and a new initial state
  vector. If no adjustments to shocks or initial state vectors are
  necessary under the policy rule, this field may be omitted.
"""
immutable AltPolicy
    rule::Function
    forecast_init::Function
end

function AltPolicy(rule)
    AltPolicy(rule, identity)
end
