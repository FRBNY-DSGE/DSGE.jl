# Alternative Policies

``` @meta
CurrentModule = DSGE
```

## Procedure

This section describes forecasting under an alternative monetary policy
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

## The `AltPolicy` Type

```@docs
DSGE.AltPolicy
```
