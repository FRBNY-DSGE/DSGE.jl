"""
```
type AltPolicy
```

Type defining an alternative policy rule.

### Fields

- `key::Symbol`: alternative policy identifier

- `eqcond::Function`: a version of `DSGE.eqcond` which computes the equilibrium
  condition matrices under the alternative policy. Like `eqcond`, it should take
  in one argument of type `AbstractModel` and return the `Γ0`, `Γ1`, `C`, `Ψ`,
  and `Π` matrices.

- `solve::Function`: a version of `DSGE.solve` which solves the model under the
  alternative policy. Like `DSGE.solve`, it should take in one argument of type
  `AbstractModel` and return the `TTT`, `RRR`, and `CCC` matrices.

- `setup::Function`

- `forecast_init::Function`: a function that initializes forecasts under the
  alternative policy rule. Specifically, it accepts a model, an `nshocks` x
  `n_forecast_periods` matrix of shocks to be applied in the forecast, and a
  vector of initial states for the forecast. It must return a new matrix of
  shocks and a new initial state vector. If no adjustments to shocks or initial
  state vectors are necessary under the policy rule, this field may be omitted.

- `color::Colorant`: color to plot this alternative policy in. Defaults to blue.

- `linestyle::Symbol`: line style for forecast plots under this alternative
  policy. See options from `Plots.jl`. Defaults to `:solid`.

"""
type AltPolicy
    key::Symbol
    eqcond::Function
    solve::Function
    setup::Function
    forecast_init::Function
    color::Colorant
    linestyle::Symbol
end

function AltPolicy(key::Symbol, eqcond_fcn::Function, solve_fcn::Function;
                   forecast_init::Function = identity,
                   setup::Function = identity,
                   color::Colorant = RGB(0., 0., 1.),
                   linestyle::Symbol = :solid)

    AltPolicy(key, eqcond_fcn, solve_fcn, setup, forecast_init, color, linestyle)
end

Base.string(a::AltPolicy) = string(a.key)