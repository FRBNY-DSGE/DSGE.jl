# [Alternative Policies](@id altpol-doc)

``` @meta
CurrentModule = DSGE
```

This section describes forecasting under an alternative monetary policy
rule. That is:

1. The state space system is split into two exogenous and unanticipated regimes,
the "historical" regime and the "forecast" regime. The historical policy rule applies
during the "historical" regime, and the alternative policy rule applies
to the "forecast' regime.
2. Filtering and smoothing is done under the historical monetary policy rule,
   i.e. the one defined in the `eqcond` method for the given model.
3. Forecasts and IRFs are computed under the alternative rule.

See [Regime-Switching Forecasts](@id regime-switch-forecast) for details on how
forecasting works when the state space system includes exogenous regime-switching.

Alternative policies can be either permanent or temporary. To use alternative policies,
the user needs to ensure that the model can be solved with regime-switching (see [Regime-Switching](@ref solveregswitch)).

## Procedure for Permanent Alternative Policies

The user defines some instance of the `AltPolicy` type (described below) and
then calls the function `setup_permanent_altpol!`. Then the function
calls made to forecast and compute means and bands remain the same as usual (see
[Forecasting](@ref) and [Computing Means and Bands](@ref means-bands)).

For example, suppose you have defined the functions `taylor93_eqcond` and
`taylor93_solve` corresponding to
[Taylor (1993)](http://www.sciencedirect.com/science/article/pii/016722319390009L)'s
proposed monetary policy rule. Then you can run:

```julia
m = AnSchorfheide()
setup_permanent_altpol!(m, AltPolicy(:taylor93, taylor93_eqcond, taylor93_solve); cond_type = :none)
forecast_one(m, :mode, :none, [:forecastobs, :forecastpseudo])
compute_meansbands(m, :mode, :none, [:forecastobs, :forecastpseudo])
```

Permanent alternative policies utilize some of the same
machinery as temporary alternative policies, but they use different algorithms
for converting the equilibrium conditions from gensys form to the
reduced form transition matrices for a state space system.
The function `setup_permanent_altpol!` performs the setup required
to interface with this machinery. The keyword argument `cond_type` is necessary
because when the alternative policy is applied depends on whether
the forecast is conditional or not. If a forecast is conditional,
it is assumed that the alternative policy does not occur until
after the conditional horizon, to maintain the idea that the alternative policy
is entirely unanticipated.

## [Procedure for Temporary Alternative Policies](@id tempaltpol-procedure)

Another counterfactual exercise is temporarily imposing a different monetary policy
rule, i.e. a temporary alternative policy, before switching
back to either the historical rule or some (permanent) alternative rule. To implement this, we utilize
exogenous regime switching in the forecast horizon.

In a rational expectations equilibrium, agents will take into account the fact that
the temporary policy is expected to terminate. A different algorithm than Chris Sims's
standard `gensys` algorithm is required, which we have implemented as `gensys2`. Note
that this `gensys2` is different from the `gensys2` Chris Sims has implemented to
calculate second-order perturbations.

To set up a temporary alternative policy, a user first needs to specify
alternative policy using the type `AltPolicy`. For instance, this code implements a
[Nominal GDP targeting policy](https://github.com/FRBNY-DSGE/DSGE.jl/blob/main/src/altpolicy/ngdp_target.jl),
and the `AltPolicy` is constructed by calling `DSGE.ngdp()`, or equivalently

```
AltPolicy(policy, DSGE.ngdp_eqcond, DSGE.ngdp_solve, forecast_init = DSGE.ngdp_forecast_init)
```

where the inputs to `AltPolicy` here are

```@docs
DSGE.ngdp_eqcond
DSGE.ngdp_replace_eq_entries
DSGE.ngdp_solve
DSGE.ngdp_forecast_init
```

Note that `ngdp_replace_eq_entries` is called by `ngdp_eqcond` but is not a direct input to `AltPolicy`.

The user also needs to complete the following steps to apply temporary alternative policies.

- Adding a regime for every period during which the alternative policy applies,
  plus one more regime for the policy which will be permanently in place after the temporary policies end.

- Adding the setting `Setting(:gensys2, true)` to indicate `gensys2` should be used. If this setting is false
  or non-existent, then alternative policies will be treated as if they are permanent. Their equilibrium
  conditions will be solved using `gensys`, which can lead to determinacy and uniqueness problems if
  the alternative policy should be temporary (e.g. a temporary ZLB).

- Adding the setting `Setting(:replace_eqcond, true)` to indicate equilibrium conditions will be replaced.

- Adding the setting `Setting(:regime_eqcond_info, info)`, where `info` should be a
  `Dict{Int, DSGE.EqcondEntry}` mapping regimes to instances of `EqcondEntry`, a type which holds
  any information needed to update equilibrium conditions to implement a given alternative policy.
  Borrowing the example of temporary NGDP targeting, the relevant `EqcondEntry` would be constructed as
  `EqcondEntry(DSGE.ngdp())`. Note that the user only needs to populate this dictionary with regimes in
  which the `eqcond` function differs from the default.

To see an example of using temporary alternative policies, see the
[example script for regime-switching](https://github.com/FRBNY-DSGE/DSGE.jl/blob/main/examples/regime_switching.jl).

## [Alternative Policy Uncertainty and Imperfect Awareness](@ref uncertainaltpol)
Click on the section header for details on how to add policy uncertainty or
imperfect credibility to alternative policies (both permanent and temporary).

## [`MultiPeriodAltPolicy`](@ref imperfect-awareness-multiperaltpol)
Click on the section header to see the primary use of the type `MultiPeriodAltPolicy`,
which extends `AltPolicy` to specify multiple regimes. In particular,
one of the fields of `MultiPeriodAltPolicy` is `regime_eqcond_info`,
which stores a dictionary that can be used to update a model's
`regime_eqcond_info` `Setting`, i.e.

```
m <= Setting(:regime_eqcond_info, multi_period_altpol.regime_eqcond_info)
```

## [Types](@id altpol-types)

```@docs
DSGE.AltPolicy
DSGE.EqcondEntry
DSGE.MultiPeriodAltPolicy
```
