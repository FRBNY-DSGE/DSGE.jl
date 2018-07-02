# Alternative Scenarios

``` @meta
CurrentModule = DSGE
```

An alternative scenario is a set of forecasted paths for a subset of observables
("targets"), along with a subset of shocks ("instruments") which get the model
to hit those paths. Each scenario has a description which tells a story about
the scenario targets and instruments -- for example, "High Spreads" or
"Persistent Consumer Optimism".

In an alternative scenario, we first solve for the values for specific shocks at
specific horizons that create the paths of the observables imposed in the
scenario. Only then do we run a regular forecast, imposing the path of shocks we
just solved for and letting these propagate through the economy. (Compare to
adding more periods of conditional data -- in this case, the model is free to
use as many shocks as it wants to explain the deviation of the conditional data
from trend in these additional periods.) This allows us to answer questions
like:

- Suppose a spread shock hits tomorrow, increasing spreads by 50 basis
  points. How large of a spread shock did this correspond to it? How does the
  economy respond? Suppose it is a discount factor shock instead.
- Suppose a discount factor shock hits for the next two quarters, causing a 100
  basis point increase in spreads for two quarters, then reverting to zero
  (increase from baseline).  How does the economy respond?

We will conceive of all alternative scenario forecasts as occuring in
"deviations from baseline forecast". That is, if the scenario consists of
spreads increasing by 50 basis points, then we want this to be independent of
the underlying model and baseline forecast, such that in our own forecast,
spreads actually are 50 basis points above our baseline forecast.


## Basic Scenarios

### The `Scenario` Type

The `Scenario` type encodes the basic information needed to forecast an
alternative scenario. For now, we'll consider the following fields of the
`Scenario` type:

- `key::Symbol`: scenario identifier, appearing in file names
- `description::String`: longer name, e.g. "High Spreads"
- `target_names::Vector{Symbol}`: names of the observables targeted in this
  scenario. For the High Spreads scenario, this would just be
  `[:obs_spread]`. These observable names should correspond to keys in
  `m.observables`
- `instrument_names::Vector{Symbol}`: names of the shocks used to hit the target
  paths, which should correspond to keys in `m.exogenous_shocks`. There must be
  at least as many instruments as there are targets. If this field is the empty
  array, then all model shocks will be used
- `targets::DataFrame`: contains the specific target values, *given in deviations
  from some baseline forecast*. No `:date` field is required, since the scenario
  is assumed to begin in the first forecasted period of the model
- `instruments::DataFrame`: initially empty `DataFrame`, which is populated
  after backing out the necessary shocks to hit the targets
- `vintage::String`: scenario vintage in `yymmdd` format, which can be different
  from the data vintage
- `n_draws::Int`: number of scenario draws. A scenario draw is one set of target
  paths. Often, multiple draws of target paths will be associated with one
  scenario name, which collectively make up a forecast distribution for the
  particular scenario. This field is usually initialized to 0 and then updated
  upon reading in the target draws

### Setting Up Input Target Paths

If you only want to forecast one scenario draw, it is sufficient to create an
instance of the `Scenario` type and hard-code in the target paths. However, when
forecasting multiple draws, it is convenient to load target draws from a file.

This file's name should be of the form
`inpath(m, "scenarios", <key>_<vintage>.jld)` and should contain the following
datasets:

- `arr::Array{Float64, 3}`: array of target values of size `ndraws` x `ntargets`
  x `horizon`
- `target_indices::OrderedDict{Symbol, Int}`: maps target names to their indices
  along the `ntargets` dimension

### Forecasting Scenarios

Forecasting scenarios is similar to running a normal full-distribution forecast,
with some exceptions:

- All draws of a particular scenario are forecasted under the modal parameters.
- Smoothed histories, shock decompositions, and IRFs are not supported.
- We zero out the entries in the ``Q`` matrix (the variance-covariance of the
  shocks ``\epsilon_t``) corresponding to shocks which are not scenario
  instruments.
- Forecasting is done in deviations from baseline. That is, let ``s^a_t`` and
  ``s^b_t`` be the state vectors under the alternative and baseline scenarios
  respectively, and define ``y^a_t`` and ``y^b_t`` analogously for observable
  vectors. Then the state space in deviations is

  ```math
  \begin{align*}
  s^a_t - s^b_t &= T(s^a_{t-1} - s^b_{t-1}) + R \epsilon_t & \mathrm{(transition)} \\
  y^a_t - y^b_t &= Z(y^a_t - y^b_t) & \mathrm{(measurement)}
  \end{align*}
  ```

The main function to run is `forecast_scenario`, which is similar in spirit to
`forecast_one`. `forecast_scenario` loads the modal parameters for the model and
calls `forecast_scenario_draw`, which does the following for each input draw:

1. Call `load_scenario_targets!` to read the `i`th scenario draw from the input
   file
2. Filter and smooth to back out the necessary shocks, treating the targeted
   paths like data
3. Forecast paths for all variables (not just the targeted observables) using
   the smoothed shocks
4. Check that the forecasted observables match the targets and return

When all draws for the scenario have been forecasted, the forecasted observables
and pseudo-observables are written to the files given by
`get_scenario_output_files`.

### Computing Means and Bands

Computing means and bands is carried out by `scenario_means_bands`, which is
likewise similar to `compute_meansbands` for regular forecasts. The default
`output_vars` passed into `scenario_means_bands` are

```
output_vars = [:forecastutobs, :forecastobs, :forecast4qobs,
               :forecastutpseudo, :forecastpseudo, :forecast4qpseudo]
```

The product `:forecastut` refers to untransformed forecasts, i.e. forecasts in
model units (in deviations from baseline).

Since these forecasts are given in deviations from baseline, different
transformations are used than in the usual case. These don't add back population
growth to per-capita variables and approximate annualizing log differences
(exponentiating and raising to the fourth power) by multiplying by four. The
mapping from usual to scenario transformation is given by
`get_scenario_transform`.


## Additional Features

### Other `Scenario` Fields

The remaining fields in the `Scenario` type are used to forecast scenarios using
additional bells and whistles:

- `shock_scaling::Float64`: after filtering and smoothing shocks, multiply them
  by `shock_scaling` before forecasting. Defaults to `1.0`.
- `draw_states::Bool`: if `true`, use the simulation smoother rather than the
  Kalman smoother to smooth shocks. This generates uncertainty around the
  target paths. Defaults to `false`.
- `altpolicy::AltPolicy`: solve for shocks under the historical policy rule,
  then switch to the alternative policy (see [Alternative Policies](@ref))
  before forecasting. Defaults to `AltPolicy(:historical, eqcond, solve)`.

### `SwitchingScenario`s

Suppose you want to simulate a set of scenario draws which all start from some
default scenario, switch to an alternative scenario (which we call the original
scenario) with some probability in each period, and then revert back to the
default scenario with some probability. This functionality is encoded in the
`SwitchingScenario` type, which has the following fields:

- `key::Symbol`: identifier, not the same as the original scenario's key
- `description::Symbol`: defaults to the original scenario's description
- `vintage::String`: defaults to the original scenario's vintage
- `original::SingleScenario`
- `default::SingleScenario`
- `probs_enter::Vector{Float64}`: gives the probability of leaving the default
  scenario and entering the original scenario in each period
- `probs_exit::Vector{Float64}`: gives the probability of leaving the original
  scenario and reverting to the default scenario in each period

`SingleScenario` is the abstract supertype of both `Scenario` and
`SwitchingScenario`.

`SwitchingScenario`s cannot be forecasted like ordinary `Scenario`s. Instead,
after simulating draws from both the original and default scenarios, you call
`simulate_switching` to use these already-forecasted draws. For each draw `i` of
the original scenario:

1. Randomly select a draw `j` of the default scenario.
2. In each period `t` (beginning in period 1), determine whether to switch to
   the original scenario with probability `probs_enter[t]`. If remaining in the
   default scenario, use the period `t` forecasted values of the `j`th default
   scenario draw.
3. Suppose the last period in the default scenario is `t0`. Then for each
   subsequent period `t0 + h`, decide whether to revert to the default scenario
   with probability `probs_exit[t0 + h]`. If remaining in the original scenario,
   use the period `h` (not `t0 + h`) forecasted values of the `i`th original
   scenario draw.
4. Let `t1` be the last period in the original scenario. Then for each remaining
   period `t1 + h`, use the period `t1 + h` forecasted values of the `j`th
   default scenario draw.

`simulate_switching` saves simulated draws in the same format as
`forecast_scenario`. Transforming and computing means and bands for
`SwitchingScenario`s using `scenario_means_bands` is the same as for regular
`Scenario`s.


## Aggregating Multiple Scenarios

Finally, we are sometimes interested in aggregating the forecast draws from
multiple scenarios. We define the `ScenarioAggregate` type, which has the
following fields:

- `key::Symbol`
- `description::String`
- `scenarios::Vector{AbstractScenario}`: vector of component scenarios, **some
  of which might be themselves `ScenarioAggregates`**
- `sample::Bool`: indicates whether to
- `proportions::Vector{Float64}`: vector of relative scenario proportions
- `total_draws::Int`: desired final number of draws
- `replace::Bool`: indicates whether to sample with replacement
- `vintage::String`

In addition to the default constructor, there are two more `ScenarioAggregate`
constructors, corresponding to the two possible values of `sample`. The key,
description, vector of component scenarios, and vintage are always specified.

- `sample = false`: No additional fields are required for this constructor, as
  the component scenario draws are kept in their original proportions. The
  `proportions` and `total_draws` fields are initialized to dummy values and are
  updated when the component draws are read in. `replace` is set to false.
- `sample = true`: Additionally specify `proportions`, `total_draws`, and
  `replace`. Draws from `scenarios[i]` are sampled into the aggregate
  distribution with probability `proportions[i]`

`SingleScenario` and `ScenarioAggregate` are both subtypes of the abstract type
`AbstractScenario`. The actual sampling and aggregating of scenario draws
happens entirely in `scenario_means_bands`, which calls functions that dispatch
on the `AbstractScenario` subtype passed in. In particular, this means that we
don't save a separate "raw" output file with the individual draws that made it
into a particular `ScenarioAggregate` -- we only save the resulting
`MeansBands`.
