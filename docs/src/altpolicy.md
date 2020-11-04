# Alternative Policies

``` @meta
CurrentModule = DSGE
```

## Procedure for Permanent Alternative Policies

This section describes forecasting under a (permanent) alternative monetary policy
rule. That is:

1. Filtering and smoothing is done under the historical monetary policy rule,
   i.e. the one defined in the `eqcond` method for the given model.
2. Before forecasting, the state space matrices are recomputed under the
   alternative policy rule.
3. Forecasts and IRFs are computed under the alternative rule.

Note that shock decompositions (and the two associated products, trends and
deterministic trends) cannot currently be computed under an alternative policy.

The user defines some instance of the `AltPolicy` type (described below) and
sets it as the value for the `:alternative_policy` `Setting`. Then the function
calls made to forecast and compute means and bands remain the same as usual (see
[Forecasting](@ref) and [Computing Means and Bands](@ref means-bands)).

For example, suppose you have defined the functions `taylor93_eqcond` and
`taylor93_solve` corresponding to
[Taylor (1993)](http://www.sciencedirect.com/science/article/pii/016722319390009L)'s
proposed monetary policy rule. Then you can run:

```julia
m = AnSchorfheide()
m <= Setting(:alternative_policy, AltPolicy(:taylor93, taylor93_eqcond, taylor93_solve))
forecast_one(m, :mode, :none, [:forecastobs, :forecastpseudo])
compute_meansbands(m, :mode, :none, [:forecastobs, :forecastpseudo])
```

## Procedure for Temporary Alternative Policies

Another counterfactual exercise is temporarily imposing a different monetary policy
rule, i.e. a temporary alternative policy. To implement this, we utilize
exogenous regime switching in the forecast horizon. See [Regime-Switching Forecasts](@id regime-switch-forecast)
for details on regime-switching.

In a rational expectations equilibrium, agents will take into account the fact that
the temporary policy is expected to terminate. A different algorithm than Chris Sims's
standard `gensys` algorithm is required, which we have implemented as `gensys_cplus`.

To set up a temporary alternative policy, a user needs to specify
the changes to the equilibrium conditions to the policy rule in `eqcond`
For instance, a
[Nominal GDP targeting policy](https://github.com/FRBNY-DSGE/DSGE.jl/blob/master/src/altpolicy/ngdp_target.jl).
uses the function `ngdp_replace_eq_defines` entries to define these changes.

```@docs
DSGE.ngdp_replace_eq_entries
```

The user also needs to complete the following steps.

- Adding a regime for every period in the forecast horizon during which the alternative policy applies,
  plus one more regime for the first regime in which the alternative policy does NOT apply.
- Adding the setting `Setting(:gensys2, true)` to indicate `gensys_cplus` should be used
- Adding the setting `Setting(:replace_eqcond, true)` to indicate the `eqcond` function will be replaced
- Adding the setting `Setting(:replace_eqcond_func_dict, replace_eqcond)`, where `replace_eqcond`
  should be a `Dict{Int, Function}` mapping regimes to alternative `eqcond` functions. Note that
  the user only needs to populate regimes in which the `eqcond` function differs from the standard one.

To see an example of using temporary alternative policies, see the
[example script for regime-switching](https://github.com/FRBNY-DSGE/DSGE.jl/blob/master/examples/regime_switching.jl).

In contrast, the permanent version of Nominal GDP targeting would be

```
AltPolicy(policy, DSGE.ngdp_eqcond, DSGE.ngdp_solve, forecast_init = DSGE.ngdp_forecast_init)
```

where

```@docs
DSGE.ngdp_eqcond
DSGE.ngdp_solve
DSGE.ngdp_forecast_init
```


## The `AltPolicy` Type

```@docs
DSGE.AltPolicy
```
