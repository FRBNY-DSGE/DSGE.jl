"""
```
mutable struct AltPolicy
```

Type defining an alternative policy rule.

### Fields

- `key::Symbol`: alternative policy identifier

- `eqcond::Function`: a version of `DSGE.eqcond` which computes the equilibrium
  condition matrices under the alternative policy. Like `DSGE.eqcond`, it should
  take in one argument of mutable struct `AbstractDSGEModel` and return the `Γ0`,
  `Γ1`, `C`, `Ψ`, and `Π` matrices.

- `solve::Function`: a version of `DSGE.solve` which solves the model under the
  alternative policy. Like `DSGE.solve`, it should take in one argument of mutable
  struct `AbstractDSGEModel` and return the `TTT`, `RRR`, and `CCC` matrices.

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
mutable struct AltPolicy
    key::Symbol
    eqcond::Function
    solve::Function
    forecast_init::Function
    color::Colorant
    linestyle::Symbol
end

function AltPolicy(key::Symbol, eqcond_fcn::Function, solve_fcn::Function;
                   forecast_init::Function = identity,
                   color::Colorant = RGB(0., 0., 1.),
                   linestyle::Symbol = :solid)

    AltPolicy(key, eqcond_fcn, solve_fcn, forecast_init, color, linestyle)
end

Base.string(a::AltPolicy) = string(a.key)

function altpolicy_replace_eq_entries(key::Symbol)
    if key == :flexible_ait
        return flexible_ait_replace_eq_entries
    elseif key == :zero_rate
        return zero_rate_replace_eq_entries
    elseif key == :smooth_ait_gdp_alt
        smooth_ait_gdp_alt_replace_eq_entries
    elseif key == :smooth_ait_gdp
        smooth_ait_gdp_replace_eq_entries
    elseif key == :ait
        ait_replace_eq_entries
    elseif key == :ngdp
        ngdp_replace_eq_entries
    elseif key == :rw_zero_rate
        rw_zero_rate_replace_eq_entries
    elseif key == :rw
        rw_replace_eq_entries
    else
        error("AltPolicy $(key) has not been added to `altpolicy_replace_eq_entries` as having a method for replacing equilibrium conditions")
    end
end

# Type to hold the entries in the regime_eqcond_info dictionary for alternative policies & regime switching
mutable struct EqcondEntry
    alternative_policy::Union{AltPolicy, Missing}
    weights::Union{Array{Float64, 1}, Missing}
end

function EqcondEntry()
    return EqcondEntry(missing, missing)
end
