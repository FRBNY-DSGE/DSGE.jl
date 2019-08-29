# Impulse Responses

```@meta
CurrentModule = DSGE
```

The forecast step allows the user to automatically compute impulse responses,
but for some purposes, a user may just want impulse responses
without having to compute any other quantities. We provide this functionality
with the `impulse_responses` function:

```@docs
DSGE.impulse_responses
```

We overload impulse_responses to cover specific use cases. For any
`AbstractDSGEModel`, we can compute the impulse responses over
a specified horizon
for all endogenous state variables, observables, and pseudo-observables
by running

```julia
m = AnSchorfheide()
system = compute_system(m)
horizon = 40
states_irf, obs_irf, pseudo_irf = impulse_response(system, horizon)
``

For an `AbstractRepModel` (a subtype of `AbstractDSGEModel` for representative
agent models), we can also grab the impulse responses by running


```julia
states_irf, obs_irf, pseudo_irf = impulse_response(m, system)
```

This use case requires the user to add a setting
under the key `:impulse_response_horizons`, which is
set by default to 40.


If a user wants to specify a subset of the exogenous shocks
and the size of those shocks, the user can run

```julia
shock_names = [:g_sh, :b_sh]
shock_values = [1.0, 1.0]
impulse_responses(m, system, horizon, shock_names, shock_values)
```

For the response of an endogenous state or observable to a specific shock,

```julia
shock_name  =  :g_sh
var_name = :obs_gdp
var_value = 0.
impulse_responses(m, system, horizon, shock_name , var_name, var_value)
```
