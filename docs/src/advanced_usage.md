# [Advanced Usage](@id advanced-usage)

```@meta
CurrentModule = DSGE
```

## Package Directory Structure

The package directory structure follows Julia module conventions. Directories in square brackets indicate future additions. *Note that this directory tree is not linked, although it appears to be.*

```@contents
Pages = ["pkg_structure.md"]
Depth = 5
```

## Working with Settings

There are many computational settings that affect how the code runs without affecting the
mathematical definition of the model. While the default settings loaded are
intended to be comprehensive rather than the minimal number of settings, users
will generally want to check that these three settings are properly chosen:

- `saveroot::String`: The root directory for model output.
- `dataroot::String`: The root directory for model input data.
- `data_vintage::String`: Data vintage, formatted `yymmdd`. By default,
  `data_vintage` is set to today's date. It is (currently) the only setting
  printed to output filenames by default.

Many functions in DSGE.jl will either require input data or create output data,
so it is important to check that the saveroot and dataroot are set as the user intends.
Setting the data vintage is also useful for reproducibility. Economic data like GDP
are frequently revised, which can pose issues for reproducing results. Setting
the data vintage allows users to guarantee the correct vintage of data is used
when generating results. By default, the data vintage is set to the current date,
so a user will need to manually set the data vintage to the desired date.

Below, we describe several important settings for package usage.

For more details on implementation and usage of settings, see [ModelConstructors.jl](https://github.com/FRBNY-DSGE/ModelConstructors.jl).

See [defaults.jl](https://github.com/FRBNY-DSGE/DSGE.jl/blob/main/src/defaults.jl) for the complete description of default settings.

#### General

- `saveroot::String`: The root directory for model output.
- `use_parallel_workers::Bool`: Use available parallel workers in computations.
- `nominal_rate_observable`: Name (as a `Symbol`) of the observable used to measure the nominal interest rate
  used to implement monetary policy.
- `monetary_policy_shock`: Name (as a `Symbol`) of the exogenous monetary policy shock in a concrete subtype of `AbstractDSGEModel`.
- `n_mon_anticipated_shocks::Int`: Number of anticipated policy shocks.
- `antshocks::Dict{Symbol, Int}`: a dictionary mapping the name of an anticipated shock
  to the number of periods of anticipation, e.g. `:b => 2` adds anticipated `b` shocks up to
  two periods ahead.
- `ant_eq_mapping::Dict{Symbol, Symbol}`: a dictionary mapping the name of an anticipated shock
  to the name of the state variable in the equation defining the shock's exogenous process, e.g. `:b => :b`
  maps an anticipated `b` shock to the equation `eq_b`.
- `ant_eq_E_mapping::Dict{Symbol, Symbol}`: a dictionary mapping the name of an anticipated shock
  to the name of the state variable in the equation defining the shock's one-period ahead expectation
  e.g. `:b => :Eb`  maps an anticipated `b` shock to the equation `eq_Eb`, where `Eb` is ``E_t[b_{t + 1}]``.
- `proportional_antshocks::Vector{Symbol}`: a vector of the names of one-period ahead anticipated shocks which are specified
  as directly proportional to the realizations of the current period's unanticipated shocks. For a shock `b`,
  the automatically generated parameter `σ_b_prop` defines the proportionality to the current period shock, e.g.
  a value of 1 indicates an anticipated shock in the next period of the same size as the current period's unanticipated shock.

#### Data and I/O

- `dataroot::String`: The root directory for model input data.
- `data_vintage::String`: Data vintage, formatted `yymmdd`. By default,
  `data_vintage` is set to today's date. It is (currently) the only setting
  printed to output filenames by default.
- `dataset_id::Int`: Dataset identifier. There should be a unique dataset ID for
  each set of observables.
- `cond_vintage::String`: Conditional data vintage, formatted `yymmdd`.
- `cond_id::Int`: Conditional dataset identifier. There should be a unique
  conditional dataset ID for each set of input, raw data mnemonics (not
  observables!).
- `cond_semi_names::Vector{Symbol}` and `cond_full_names::Vector{Symbol}`: names
  of observables for which we want to use semi- and full conditional data. All
  other observables are `NaN`ed out in the conditional data periods.
- `population_mnemonic::Nullable{Symbol}`: population series mnemonic in form
  `Nullable(:<mnemonic>__<source>)` (for example, `Nullable(:CNP16OV__FRED)`),
  or `Nullable{Symbol}()` if the model doesn't use population data

#### Dates

- `date_presample_start::Date`: Start date of pre-sample.
- `date_mainsample_start::Date`: Start date of main sample.
- `date_zlb_start::Date`: Start date of zero lower bound regime.
- `date_zlb_end::Date`: End date of zero lower bound regime.
- `date_forecast_start::Date`: Start date of forecast period (or the period
  after the last period for which we have GDP data).
- `date_forecast_end::Date`: End date of forecast, i.e. how far into the future to forecast.
- `date_conditional_end::Date`: Last date for which we have conditional
  data. This is typically the same as `date_forecast_start` when we condition on
  nowcasts and current quarter financial data.


#### Estimation

##### Metropolis-Hastings Settings
- `reoptimize::Bool`: Whether to reoptimize the posterior mode. If `true` (the
    default), `estimate` begins reoptimizing from the model object's parameter
    vector.  See [Optimizing or Reoptimizing](@ref estimation-reoptimizing) for
    more details.
- `calculate_hessian::Bool`: Whether to compute the Hessian. If `true` (the
    default), `estimate` calculates the Hessian at the posterior mode.
- `n_mh_simulations::Int`: Number of draws from the posterior distribution per
  block.
- `n_mh_blocks::Int`: Number of blocks to run Metropolis-Hastings.
- `n_mh_burn::Int`: Number of blocks to discard as burn-in for Metropolis-Hastings.
- `mh_thin::Int`: Metropolis-Hastings thinning step.
- `parallel::Bool`: Flag for running algorithm in parallel.
- `n_parts::Int`: Number of particles.
- `n_blocks::Int`: Number of parameter blocks in mutation step.
- `n_mh_steps::Int`: Number of Metropolis Hastings steps to attempt during the mutation step.

##### Sequential Monte Carlo Settings
- `λ::S`: The 'bending coefficient' λ in Φ(n) = (n/N(Φ))^λ
- `n_Φ::Int`: Number of stages in the tempering schedule.
- `resampling_method::Symbol`: Which resampling method to use.
    - `:systematic`: Will use sytematic resampling.
    - `:multinomial`: Will use multinomial resampling.
    - `:polyalgo`: Samples using a polyalgorithm.
- `threshold_ratio::S`: Threshold s.t. particles will be resampled when the population
    drops below threshold * N
- `c::S`: Scaling factor for covariance of the particles. Controls size of steps in mutation step.
- `α::S`: The mixture proportion for the mutation step's proposal distribution.
- `target::S`: The initial target acceptance rate for new particles during mutation.

- `use_chand_recursion::Bool`: Flag for using Chandrasekhar Recursions in Kalman filter.
- `use_fixed_schedule::Bool`: Flag for whether or not to use a fixed tempering (ϕ) schedule.
- `tempering_target::S`: Coefficient of the sample size metric to be targeted when solving
    for an endogenous ϕ.
- `old_data::Matrix{S}`:
- `old_cloud::ParticleCloud`:
- `old_vintage::String`: String for vintage date of old data
- `smc_iteration::Int`: The iteration index for the number of times SMC has been run on the
     same data vintage. Primarily for numerical accuracy/testing purposes.

- `run_test::Bool`: Flag for when testing accuracy of program
- `filestring_addl::Vector{String}`: Additional file string extension for loading old cloud.
- `continue_intermediate::Bool`: Flag to indicate whether one is continuing SMC from an
    intermediate stage/
- `intermediate_stage_start::Int`: Intermediate stage at which one wishes to begin the estimation.
- `save_intermediate::Bool`: Flag for whether one wants to save intermediate Cloud objects
- `intermediate_stage_increment::Int`: Save Clouds at every increment
    (1 = each stage, 10 = every 10th stage, etc.)


#### Forecasting

- `forecast_jstep::Int`: Forecast thinning step.
- `forecast_block_size::Int`: Number of draws in each forecast block *before*
  thinning by `forecast_jstep`.
- `forecast_input_file_overrides::Dict{Symbol, String}`: Maps `input_type`(s) to
  the file name containing input draws for that type of forecast. See
  [Forecasting](@ref).
- `forecast_horizons::Int`: Number of periods to forecast.
- `impulse_response_horizons::Int`: Number of periods for which to calculate IRFs.
- `n_periods_no_shocks::Int`: Number of periods for which no shocks are drawn (e.g.
  a full-distribution forecast draws shocks, but if `n_periods_no_shocks = 3`, then for
  3 periods in the forecast horizon, no shocks will be drawn)

#### Alternative Policy

- `alternative_policy::AltPolicy`: See [Alternative Policies](@ref).

### Accessing Settings

The function `get_setting(m::AbstractModel, s::Symbol)` returns the value of the setting `s`
in `m.settings`. Some settings also have explicit getter methods that take only the model
object `m` as an argument. Note that not all are exported.

### Overwriting Default Settings

To overwrite default settings added during model construction, a user must
create a `Dict{Symbol, Setting}` and pass that into the model constructor as the
keyword argument `custom_settings`. If the `print`, `code`, and `description`
fields of the new `Setting` object are not provided, the fields of the existing
setting will be maintained. If new values for `print`, `code`, and `description`
are specified, and if these new values are distinct from the defaults for those
fields, the fields of the existing setting will be updated.

For example, overwriting `use_parallel_workers` should look like this:
```julia
custom_settings = Dict{Symbol, Setting}(
    :use_parallel_workers => Setting(:use_parallel_workers, true))
m = Model990(custom_settings = custom_settings)
```

Or like this:
```julia
m = Model990()
m <= Setting(:use_parallel_workers, true)
```
Note that using this second method will not work for all settings, e.g. `n_anticipated_shocks`
is a setting that must be passed into the model during construction, as in the first example.

By default, passing in `custom_settings` overwrites the entries in the model
object's `settings` field. However, with the additional keyword argument
`testing = true`, it will overwrite the entries in `test_settings`:

```julia
m = Model990(custom_settings = custom_settings, testing = true)
```

## [Accelerating Computation of Regime-Switching System](@id accel-regime-switch-statespace-comp)
Regime-switching state space systems take more time to compute, which can
severely slow down estimation and forecasting time. The interface for
computing regime-switching systems is written to be easy to use and generic,
but, as a result, its default behavior ignores information that could
be used to accelerate computation time. We provide some settings
that allow the user to specify such information about the state space system.

- `perfect_credibility_identical_transitions::Dict{Int, Int}`: different regimes
  may have the same transition equations (`TTT`, `RRR`, and `CCC` matrices,
  also see [Solving the Model](@ef solving-dsge-doc)).
  This setting tells the code to use another regime's transition equation rather than
  recalculate the equation. The keys of this `Dict` are regime numbers,
  and the values specify the regime to which the keys' regimes are identical.
  For example, if the setting was `Dict(2 => 1, 3 => 1)`, then we are saying
  that regimes 2 and 3 have the same transition equations as regime 1.

- `identical_eqcond_regimes::Dict{Int, Int}`: different regimes
  may have the same equilibrium conditions (see [Solving the Model](@ef solving-dsge-doc)).
  This setting tells the code to copy another regime's equilibrium conditions
  rather than recompute the gensys matrices. The keys of this `Dict` are regime numbers,
  and the values specify the regime to which the keys' regimes are identical.
  For example, if the setting was `Dict(2 => 1, 3 => 1)`, then we are saying
  that regimes 2 and 3 have the same equilibrium conditions as regime 1.

- `empty_measurement_equation::Vector{Bool}`: when using time-varying
  information sets and forward-looking observables, you may need to calculate
  the transition equations beyond the last period of available data.
  By default, `compute_system` will also compute the measurement equation
  for these regimes in the future, which is unnecessary if you are trying to
  estimate the model. This setting specifies which regimes
  can have an empty measurement equation (set to be a `Measurement` type
  with undefined matrices for its fields). A `false` element in the vector
  means that the measurement equation is nonempty, while a `true` element
  means an empty measurement equation. The length of the vector
  should be the same length as the number of regimes, and the indices
  of the vector corresponding to the regime number.

- `empty_pseudo_measurement_equation::Vector{Bool}`: same as `empty_measurement_equation`
  but for the pseudo-measurement equation. For estimations, you can omit all
  the pseudo-measurement equations since they are unnecessary.

## [Regime-Switching Forecasts](@id regime-switch-forecast)

Forecasts can involve state-space systems with exogenous and unanticipated regime-switching
in the history periods and forecast horizon. Anticipated temporary alternative policies
can also occur in both the history and the forecast horizon. Historical
regime switching may occur to reflect structural breaks or
to allow a DSGE to handle special circumstances, such as the COVID-19 pandemic.
Regime switches in the forecast horizon may occur because agents
expect a ZLB until some date in the future. In a rational expectations equilibrium,
agents will behave differently if they know a forecasted policy is temporary
rather than permanent. Using exogenous regime-switching along with a modified
`gensys` solution algorithm is one way of implementing this expectation.
See [Regime-Switching](@ref solveregswitch) for more details on the solution algorithm.

In this section, we will go over the interface for running
regime-switching forecasts and discuss some details of the implementation.
It is useful to also look at the posted
[example script](https://github.com/FRBNY-DSGE/DSGE.jl/blob/main/examples/regime_switching.jl) for regime-switching.
To understand how to implement your own regime-switching model, we recommend
examining the implementation of [regime-switching
equilibrium conditions](https://github.com/FRBNY-DSGE/DSGE.jl/blob/main/src/models/representative/m1002/eqcond.jl)
for `Model1002` and how it is integrated with our [solvers](https://github.com/FRBNY-DSGE/DSGE.jl/blob/main/src/solve/solve.jl).
For a guide to running permanent and/or temporary alternative policies, please see [Alternative Policies](@ref).

### Preparing a Model's Settings for Regime-Switching

Suppose we wanted to run a regime-switching forecast, where
where the regimes are 1959:Q3-1989:Q4, 1990:Q1-2019:Q3, and 2019:Q4 to the
end of the forecast horizon. The following lines are required:

```
m <= Setting(:regime_switching, true)
m <= Setting(:regime_dates, Dict{Int, Date}(1 => Date(1959, 9, 30), 2 => Date(1990, 3, 31),
                                            3 => Date(2019, 12, 31))
```

The first setting turns on regime switching. Internally, functions like `forecast_one`
will decide whether to use regime switching or not depending
on whether `get_setting(m, :regime_switching)` is true and whether
there are actually multiple regimes specified by `:regime_dates`.
The second setting is a `Dict` mapping
the regime number to the first date (inclusive) of that regime.

Before running a forecast, we must also run

```
setup_regime_switching_inds!(m; cond_type = cond_type)
```

which will automatically compute the (required) settings

- `:reg_forecast_start`: Regime of the first forecast start period
- `:reg_post_conditional_end`: Regime of the period after the last conditional forecast period
- `:n_regimes`: Number of total regimes. If this is 1, then regime switching will not occur.
- `:n_hist_regimes`: Number of regimes in the history
- `:n_fcast_regimes`: Number of regimes in the forecast horizon (including the conditional forecast)
- `:n_cond_regimes`: Number of regimes in the conditional forecast

These settings will generally depend on whether the forecast is conditional or not,
so the user needs to pass in `cond_type = :full` or `cond_type = :semi` to `setup_regime_switching_inds!`
if the user wants a forecast with correct regime-switching.

Finally, to run a full-distribution forecast with regime-switching
using `forecast_one` or `usual_model_forecast`, it is necessary to manually
construct the matrix of parameter draws and pass them as an input with the keyword `params`.
Currently, we have not fully implemented loading parameters from a saved estimation file.
For an example about how to do this, see this
[example script](https://github.com/FRBNY-DSGE/DSGE.jl/blob/main/examples/regime_switching.jl).

### [Time-Varying Information Sets](@id tvis)
In many applications with regime-switching, changes in information sets may occur.
For example, say in 1959:Q1 - 2007:Q3, people expected the Federal Reserve to always
use a Taylor-style monetary policy rule. But in 2007:Q4, people realize that the Federal
Reserve will implement a zero lower bound for `N` periods before switching back to the
Taylor-style rule. The measurement equation used in 2007:Q4 and subsequent periods need
to account for this change in the information set. In particular, quantities like
anticipated nominal rates and 10Y inflation rates involve forecasting the expected
state of the economy. If the transition matrices (e.g. `TTT`) are time-varying, and
agents' information set includes knowledge that the matrices are time-varying, then
the measurement equation should account for it. Explicitly, assume the state space
evolves according to
```math
\begin{aligned}
s_t = T_t s_{t - 1} + R_t \varepsilon_t + C_t.
\end{aligned}
```

If in period `t`, the measurement equation includes the anticipated nominal rate
in ``k`` periods, and agents know that the transition equations are time-varying over some horizon ``H``,
then we need to calculate that expectation taking into account that agents know about
the time variation in the matrices, e.g. ``T_{t + 1}, \dots, T_{t + H}``. This approach
allows for varying degrees of myopia, e.g. ``H = 0`` implies that agents do not know about
any time variation while ``H < k`` captures the case that agents only know about the
time variation up to a certain horizon forward.

To help the user write the correct measurement equation with time-varying transition equations,
we have implemented the following two functions:

```@docs
DSGE.k_periods_ahead_expectations
DSGE.k_periods_ahead_expected_sums
```

See the [measurement equation for Model 1002](https://github.com/FRBNY-DSGE/DSGE.jl/blob/main/src/model/representative/m1002/measurement.jl)
for an example of how these functions are used.

To accelerate the computation time for these functions, we have also implemented types
that create memos of the products and powers of the ``\{T_{t + k}\}_{t = 1}^k`` matrices
which are needed. See

```@docs
DSGE.ForwardExpectationsMemo
DSGE.ForwardMultipleExpectationsMemo
```

The first type is mainly an "under the hood" type for `k_periods_ahead_expectations`.
The second type is a wrapper type that constructs all the memos needed
to implement forward expectations of levels and sums in an efficient manner.
We have automated the construction of memos with the Boolean `Setting`s
`use_forward_expectations_memo` and `use_forward_expected_sum_memo`. The first `Setting`
indicates that a memo type will be used for calls to `k_periods_ahead_expectations`
and the second indicates a memo type will be used for calls to `k_periods_ahead_expected_sums`.
It is assumed by default that the last matrix in the sequence `TTTs` in
a `RegimeSwitchingSystem` is the first period in which a `TTT` matrix permanently applies (hence
we may assume that in all future periods the `TTT` is the same), but
if this is not the case, then the user needs to specify the correct regime
with the `Setting` `memo_permanent_policy_regime::Int`.

For details on how we implement a state space system with time-varying information sets,
see [The `TimeVaryingInformationSetSystem` Type](@ref tvistype). For guidance on how to use this type,
e.g. calculating forecats, see this [example script](https://github.com/FRBNY-DSGE/DSGE.jl/blob/main/examples/tvis_system.jl).

### Available Types of Regime Switching

There are three cases involving regime switching that are implemented in DSGE.jl

- Exogenous and unanticipated regime switching (e.g. unanticipated regime-switching parameters)
- Alternative policies (temporary and permanent)
- Time-varying information sets

To implement regime-switching parameters or use temporary alternative policies, see this
[example script](https://github.com/FRBNY-DSGE/DSGE.jl/blob/main/examples/regime_switching.jl) for regime-switching forecasts.
This [documentation on temporary alternative policies](@ref tempaltpol-procedure) will also be helpful.
For further details on regime-switching parameters, see
the [documentation for ModelConstructors.jl](https://github.com/FRBNY-DSGE/ModelConstructors.jl).
To implement time-varying information sets, see
this [example script](https://github.com/FRBNY-DSGE/DSGE.jl/blob/main/examples/tvis_system.jl).

Once the regime-switching settings are properly created, the syntax for running a forecast
is the same as when there is no regime-switching. See [Forecasting](@ref).

### Handling of the Zero Lower Bound (ZLB)

The New York Fed DSGE model can handle the ZLB in two ways.

In the first way, the New York Fed DSGE model treats the ZLB as a temporary alternative policy
over a pre-specified horizon. In the second way, the New York Fed DSGE model treats the ZLB as a
separate regime in which anticipated monetary policy shocks become "alive" and have positive
standard deviations. However, this second form of the ZLB is not implemented
as a separate regime. The reason is the only difference in the "pre-ZLB" and "post-ZLB"
regimes is whether or not anticipated monetary policy shocks are non-zero. For an example, see
the [smoothing code](https://github.com/FRBNY-DSGE/DSGE.jl/blob/main/src/forecast/smooth.jl)
as well as the auxiliary functions `zlb_regime_matrices` and `zlb_regime_indices`
in this [file](https://github.com/FRBNY-DSGE/DSGE.jl/blob/main/src/estimate/kalman.jl).

This approach saves computational time. Rather than creating redundant matrices,
we directly zero/un-zero the appropriate entries in the pre- and post-ZLB `QQ` matrices.
This approach also economizes on unnecessary switching,
For instance, during the calculation of shock decompositions and trends, it is unnecessary
to distinguish between the pre- and post-ZLB regimes.

## [Alternative Policy Uncertainty and Imperfect Awareness](@id uncertainaltpol)
The standard alternative policy code assumes that people completely believe the change in policy.
However, in many cases, the more realistic modeling choice is assuming some uncertainty
or imperfect awareness/credibility about the policy change. This approach can also partially address the concern
that expectations have counterfactually strong effects in standard DSGEs (e.g. the forward guidance puzzle).

### Theory
We model imperfect awareness by assuming there are ``n`` possible alternative policies
that may occur tomorrow and ``n`` probability weights assigned to each policy. Further, it is believed
that the alternative policy which occurs tomorrow will be permanent. One of the policies is
the alternative policy which is actually implemented. The function
[`gensys_uncertain_altpol`](https://github.com/FRBNY-DSGE/DSGE.jl/tree/main/src/solve/gensys_uncertain_altpol.jl)
calculates the state space transition equation implied by these beliefs.
A typical application is assuming that with probability ``p`` some alternative policy occurs tomorrow
and with probability ``1-p`` the historical policy occurs tomorrow.

Imperfect awareness can occur in multiple periods and feature time-varying credibility by assuming myopia.
For example, say agents in period ``t`` believe the central bank will implement AIT in ``t + 1`` and all
subsequent periods with probability ``p_t`` and the historical rule otherwise. After period ``t + 1`` occurs and
the central bank actually implements AIT, agents again believe that in period ``t + 2`` and all subsequent periods,
the central bank will implement AIT with probability ``p_{t + 1}`` and the historical rule otherwise.

Imperfect awareness is robust to temporary alternative policies but requires the algorithm to account
for time variation in the transition equation. In particular, we need to first compute
the entire sequence of transition equations under the temporary alternative policy with
perfect credibility. Once this sequence is available, we can then treat
the temporary alternative policy as the alternative policy which is actually implemented
and apply the same calculations described in the previous paragraph to each period of the temporary alternative policy.

### Implementation

To apply imperfect awareness, the user needs to specify the possible alternative policies
and the probability weights on these policies.

The alternative policy which is actually implemented should be
added to the `:regime_eqcond_info` dictionary as an `EqcondEntry`, e.g.

```
get_setting(m, :regime_eqcond_info)[2] = DSGE.EqcondEntry(DSGE.ngdp())
```

The other alternative policies that agents believe may occur are then added as follows:

```
m <= Setting(:alternative_policies, [altpolicy1, altpolicy2]) # altpolicy1 and altpolicy2 are AltPolicy instances
```

The user specifies the probability weights when creating the `EqcondEntry` instance
for the `:regime_eqcond_info` dictionary, e.g.

```
DSGE.EqcondEntry(DSGE.ngdp(), [p_t, 1 - p_t])
```

This approach permits time-variation in the probability weight because the user can use different `p_t`
for each regime, e.g.

```
get_setting(m, :regime_eqcond_info)[2] = DSGE.EqcondEntry(DSGE.ngdp(), [.5, .5])
get_setting(m, :regime_eqcond_info)[3] = DSGE.EqcondEntry(DSGE.ngdp(), [1., 0.])
```

Finally, before solving for the state space system or running forecasts, the user needs to add the line

```
m <= Seting(:uncertain_altpolicy, true)
```

To use imperfect awareness with a temporary altpolicy (eg. ZLB), the user needs to also add the following lines to
the model's setup:

```
m <= Setting(:uncertain_temporary_altpolicy, true)
m <= Setting(:temporary_altpolicy_length, n_zlb_regs)
```

The first line tells that a temporary altpolicy with imperfect awareness should apply.
The second line indicates the number of regimes for which the temporary altpolicy occurs. If the second line is not specified,
then it is assumed that all regimes in `get_setting(m, :regime_eqcond_info)` except the last one
are temporary altpolicy regimes. This assumption can be wrong, for example, if credibility changes after the temporary altpolicy ends.

For further guidance on adding imperfect awareness, please see the script
[uncertain_altpolicy_zlb.jl](https://github.com/FRBNY-DSGE/DSGE.jl/tree/main/examples/uncertain_altpolicy_zlb.jl).

### Forward-Looking Variables in the Measurement and Pseudo-Measurement Equations
The measurement and pseudo-measurement equations often include "forward-looking" observables, such as
the anticipated nominal interest rate and the expected average inflation rate over the next ten years.
The measurement equations for these observables are therefore affected when imperfect awareness is assumed.
Say ``ZZ_1`` and ``ZZ_2`` are the measurement equation matrices mapping states to observables under two different
monetary policy rules which may occur and that all the observables are forward-looking.
For simplicity, additionally assume that the associated ``DD_1`` and ``DD_2`` are both zero.
Because agents at the end of period ``t`` believe that either policy 1 or policy 2 occurs permanently in period ``t + 1``,
the measurement equation agents use to map states to data is just the weighted average
of the measurement matrices, i.e.

```math
\begin{aligned}
ZZ = p ZZ_1 + (1 - p) ZZ_2
\end{aligned}
```

The reason is that, conditional on alternative policy 1 occurring, the observables in ``t + 1`` should be

```math
\begin{aligned}
\mathbb{E}_t[y_{t + 1} \mid \text{policy 1}] & = ZZ_1 \mathbb{E}_t[s_{t + 1} \mid \text{policy 1}].
\end{aligned}
```

Similarly, if policy 2 occurs, then

```math
\begin{aligned}
\mathbb{E}_t[y_{t + 1} \mid \text{policy 2}] & = ZZ_2 \mathbb{E}_t[s_{t + 1} \mid \text{policy 2}].
\end{aligned}
```

The law of iterated expectations gives us the desired result.

The user does not need to worry about coding their measurement equations to account for this,
as long as the measurement equation will properly compute ``ZZ_i``,
given the policies specified in the settings `:regime_eqcond_info` and `:alternative_policies`.
DSGE.jl will handle the calculation of the convex combinations under the hood. The only setting
which users are advised to add is one that indicates which rows of ``ZZ``
are associated with forward-looking observables, e.g.

```
m <= Setting(:forward_looking_observables, [:obs_longinflation, :obs_nominalrate1])
m <= Setting(:forward_looking_pseudo_observables, [:Expected10YearNaturalRate])
```

If such a setting exists, then DSGE.jl will only calculate the weighted average for the
rows associated with these observables/pseudo-observables. Otherwise, we
compute the weighted average of the different measurement matrices. This latter approach
will always work, but it comes at the cost of unnecessary operations.

### [Imperfect Awareness with Temporary Policies as Alternative Policies]
The previous documentation generally assumes that the
alternative policies which people believe may occur are one-regime and permanent
policies. However, it is possible that agents are imperfectly aware
over alternative policies that involve temporary policies and thus require the use of `gensys2`.
The only difference the user needs to do is use a `MultiPeriodAltPolicy` type
rather than an `AltPolicy` type when populating the `Setting` `alternative_policies`.
See [Types](@ref altpol-types) for documentation on the fields of a `MultiPeriodAltPolicy`.
As an example, the code snippet below implements a temporary ZLB as the alternative policy,
assuming the existence of a regime-switching model instance `m`.

```
# Alternative Policy 1: default/historical rule
altpol1 = default_policy()

# Alernative Policy 2: ZLB starting in regime 3 and ending in regime 5, and flexible AIT starting in regime 6
new_reg_eqcond_info = Dict(3 => EqcondEntry(zlb_rule(), reg3_weights), # reg3_weights specifies whatever weights
                           4 => EqcondEntry(zlb_rule(), reg4_weights), # the user wants in regime 3, etc.
                           5 => EqcondEntry(zlb_rule(), reg5_weights),
                           6 => EqcondEntry(flexible_ait(), reg6_weights))
new_infoset         = [1:1, 2:2, [i:6 for i in 3:6]..., [i:i for i in 7:get_setting(m, :n_regimes)]...]

delete!(m.settings, :alternative_policies) # if :alternative_policies already exists, then a type error may occur

altpol2 = MultiPeriodAltPolicy(:temporary_zlb, get_setting(m, :n_regimes),
                               reg_eqcond_info, gensys2 = true,
                               temporary_altpolicy_names, [:zlb_rule],
                               temporary_altpolicy_length = 3,
                               infoset = new_infoset)

# Both AltPolicy and MultiPeriodAltPolicy are subtypes of AbstractAltPolicy
m <= Setting(:alternative_policies, DSGE.AbstractAltPolicy[altpol1, altpol2])
```

## Automatically Generating Anticipated Shocks

We have implemented some functionality for automatically adding anticipated shocks
for `Model1002`. To add these shocks, the user must pass custom settings
into the constructor using the `custom_settings` keyword. The available settings
for defining these shocks are:

- `antshocks::Dict{Symbol, Int}`: a dictionary mapping the name of an anticipated shock
  to the number of periods of anticipation, e.g. `:b => 2` adds anticipated `b` shocks up to
  two periods ahead.
- `ant_eq_mapping::Dict{Symbol, Symbol}`: a dictionary mapping the name of an anticipated shock
  to the name of the state variable in the equation defining the shock's exogenous process, e.g. `:b => :b`
  maps an anticipated `b` shock to the equation `eq_b`.
- `ant_eq_E_mapping::Dict{Symbol, Symbol}`: a dictionary mapping the name of an anticipated shock
  to the name of the state variable in the equation defining the shock's one-period ahead expectation
  e.g. `:b => :Eb`  maps an anticipated `b` shock to the equation `eq_Eb`, where `Eb` is ``E_t[b_{t + 1}]``.
- `proportional_antshocks::Vector{Symbol}`: a vector of the names of one-period ahead anticipated shocks which are specified
  as directly proportional to the realizations of the current period's unanticipated shocks. For a shock `b`,
  the automatically generated parameter `σ_b_prop` defines the proportionality to the current period shock, e.g.
  a value of 1 indicates an anticipated shock in the next period of the same size as the current period's unanticipated shock.

As an example, the following code
creates an instance of `Model1002` with anticipation of `b` shocks up to two periods ahead.

```
custom_settings = Dict{Symbol, Setting}(:antshocks => Setting(:antshocks, Dict{Symbol, Int}(:b => 2)),
                :ant_eq_mapping => Setting(:ant_eq_mapping, Dict{Symbol, Symbol}(:b => :b)))
m = Model1002("ss10"; custom_settings = custom_settings)
```

## [Automatic Endogenous ZLB Enforcement as Temporary Rule](@id auto-endo-zlb)
The user can enforce the ZLB during the forecast horizon in two ways. The default approach
uses unanticipated monetary policy shocks. Instead, the user can also use
the temporary ZLB machinery to enforce the ZLB. This enforcement is endogenous in the sense
that, conditional on a draw of shocks, we want to figure out the required length of a
temporary ZLB that will deliver non-negative interest rates throughout the horizon.

The enforcement is automated by trading off two objectives. First, we want the
length of the ZLB to be minimal so that the ZLB is not unnecessarily accommodative, unless
it is specifically desired for the ZLB to extend to at least some date.
Second, we want to maintain a reasonable computational time. Finding a minimal ZLB length
when there are multiple disconnected periods of negative interest rates would be prohibitively
expensive because expecting more periods of temporary ZLB in the future affects agents' expectations
today, and changing the number of periods of temporary ZLB in the past affects
the future evolution of states.

Instead, we endogenously enforce the ZLB only for
the first connected sequence of periods with negative interest rates
and use unanticipated monetary policy shocks for future sequences of periods with negative rates. Our algorithm proceeds as follows.

1. Forecast without any periods of temporary ZLBs (unless a minimum length is specified) and find the first connected
   sequence of periods with negative interest rates.
2. Guess a sequence of temporary ZLB regimes that cover this first sequence of periods with negative interest rates.
3. If the forecast under the temporary policy from 2 successfully enforces the ZLB over that first sequence
   and does not introduce negative rates after liftoff from the ZLB,
   then test whether shorter ZLBs will also enforce it. Otherwise, extend the sequence of temporary ZLBs
   using the same approch as 2.
4. Once we have successfully found a minimum length that guarantees non-negative rates for the first sequence of periods,
   re-run the forecast using unanticipated monetary policy shocks to enforce any other sequences of periods with negative rates.

Note that sometimes extending the sequence of temporary ZLB regimes will cause two disjoint sequences of periods with negative rates
to become contiguous, in which case we treat the two disjoint sequences as one connected sequences therafter.

To use this method, the user runs a forecast as follows
```
# (optional) maximum permitted length for temporary ZLB regimes
m <= Setting(:max_temporary_altpol_length, max_zlb_length)

# (optional) minimum permitted length for temporary ZLB regimes and the
# ZLB regimes are assumed to start in the first period of the forecast
# or the first period after the conditional horizon if a conditional forecast is run.
m <= Setting(:min_temporary_altpol_length, min_zlb_length)

forecast_one(m, input_type, cond_type, output_vars; rerun_smoother = true,
             zlb_method = :temporary_altpolicy,
             set_regime_vals_altpolicy = my_set_regime_vals_altpolicy_fnct,
             set_info_sets_altpolicy = my_set_info_sets_altpolicy_fnct,
             update_regime_eqcond_info! = my_update_regime_econd_info_fnct!,
             nan_endozlb_failures = false)
```

The keyword arguments are briefly described below. For more details,
see the docstring for `forecast_one`.

- `rerun_smoother::Bool`: needs to be true if the current sequence of temporary ZLB regimes start
  during the history or conditional horizon because changing the length
  of the temporary ZLB affects the smoothed estimate of the state at the start of the forecast.

- `zlb_method::Symbol`: set to `:temporary_altpolicy` to enforce the ZLB as a temporary policy.
  Otherwise, unanticipated monetary policy shocks will be used.

- `set_regime_vals_altpolicy::Function`: if there are regime-switching parameters,
  this function is needed to figure out what parameters should be assigned
  to the new model regimes added when extending the temporary ZLB length.

- `set_info_sets_altpolicy::Function`: if the `Setting` `tvis_information_set` is used,
  then we need to specify how to update `tvis_information_set` as new model regimes are added
  to extend the temporary ZLB length.

- `update_regime_eqcond_info!::Function`: specifies how to update `regime_eqcond_info`
  to include more or fewer regimes of temporary ZLB.

- `nan_endozlb_failures::Bool`: sometimes the ZLB cannot be enforced because rates
  are negative even when the ZLB extends throughout the entire forecast horizon or because
  the max ZLB length is reached. By default, we enforce the remainder of the forecast horizon
  with unanticipated monetary policy shocks. If this kwarg is true, we return
  `NaN`s rather than use unanticipated shocks.


TODO: add example script

## [Editing or Extending a Model](@id editing-extending-model)

Users may want to extend or edit `Model990` in a number of different ways.  The most common
changes are listed below, in decreasing order of complexity:

1. Add new parameters
2. Modify equilibrium conditions or measurement equations
3. Change the values of various parameter fields (i.e. initial `value`, `prior`,
   `transform`, etc.)
4. Change the values of various computational settings (i.e. `reoptimize`,
   `n_mh_blocks`)

Points 1 and 2 often go together (adding a new parameter guarantees a change in equilibrium
conditions), and are such fundamental changes that they increment the model specification
number and require the definition of a new subtype of `AbstractModel` (for instance,
`Model991`).  See [Model specification](@ref model-specification-mspec) for more details.

Any changes to the initialization of preexisting parameters are defined as a new model
*sub-specification*, or *subspec*. While less significant than a change to the model's
equilibrium conditions, changing the values of some parameter fields (especially priors) can
have economic significance over and above settings we use for computational purposes.
**Parameter definitions should not be modified in the model object's constructor.** First,
incrementing the model's sub-specification number when parameters are changed improves
model-level (as opposed to code-level) version control. Second, it avoids potential output
filename collisions, preventing the user from overwriting output from previous estimations
with the original parameters. The protocol for defining new sub-specifications is described
in [Model sub-specifications](@ref model-sub-specifications-msubspec).

### [Model specification (`m.spec`)](@id model-specification-mspec)

A particular model, which corresponds to a subtype of `AbstractModel`, is defined as a set
of parameters, equilibrium conditions (defined by the `eqcond` function) and measurement
equations (defined by the `measurement` function).  Therefore, the addition of new
parameters, states, or observables, or any changes to the equilibrium conditions or
measurement equations necessitate the creation of a new subtype of `AbstractModel.`

To create a new model object, we recommend doing the following:

1. Duplicate the `m990` directory within the [models](https://github.com/FRBNY-DSGE/DSGE.jl/tree/main/src/models) directory. Name the new
   directory `mXXX.jl`, where `XXX` is your chosen model specification number or string.
   Rename `m990.jl` in this directory to `mXXX.jl`.
2. In the `mXXX/` directory, change all references to `Model990` to `ModelXXX`.
3. Edit the `m990.jl`, `eqcond.jl`, and `measurement.jl` files as you see fit.  If adding
   new states, equilibrium conditions, shocks, or observables, be sure to add them to the
   appropriate list in `init_model_indices`.
4. Open the module file, `src/DSGE.jl`. Add `ModelXXX` to the list of functions to export,
   and include each of the files in `src/model/mXXX`.

### [Model sub-specifications (`m.subspec`)](@id model-sub-specifications-msubspec)

`Model990` sub-specifications are initialized by overwriting initial parameter definitions
before the model object is fully constructed. This happens via a call to `init_subspec` in
the `Model990` constructor. (Clearly, an identical protocol should be followed for new model
types as well.)

To create a new sub-specification (e.g., subspec 1) of `Model990`, edit the file
`src/models/subspecs.jl` as follows (note that this example is not actually
sub-specification `1` of `Model990`. In the source code, our sub-specification `5` is
provided as additional example.):

**Step 1.** Define a new function, `ss1`, that takes an object of type `Model990` (not
   `AbstractModel`!) as an argument. In this function, construct new parameter objects and
   overwrite existing model parameters using the `<=` syntax. For example,

```julia
function ss1(m::Model990)
    m <= parameter(:ι_w, 0.000, (0.0, .9999), (0.0,0.9999), DSGE.Untransformed(), Normal(0.0,1.0), fixed=false,
                   description="ι_w: Some description.",
                   tex_label="\\iota_w")
    m <= parameter(:ι_p, 0.0, fixed=true,
                   description= "ι_p: Some description"
                   tex_label="\\iota_p")
end
```

**Step 2.** Add an `elseif` condition to `init_subspec`:

```julia
    ...
    elseif subspec(m) == "ss1"
        return ss1(m)
    ...
```
To construct an instance of `Model990`, `ss1`, call the constructor
for `Model990` with `ss1` as an argument. For example,

```julia
m = Model990("ss1")
```

## Additional Tips
* The file `abstractdsgemodel.jl` defines numerous auxiliary functions, which allow the
  user to more easily call standard settings or count the number of dimensions for
  important variables. For example, `data_vintage(m)` returns the vintage of the data
  specified by the model object `m`. Additionally see `abstractmodel.jl` in
  [ModelConstructors.jl](https://github.com/FRBNY-DSGE/ModelConstructors.jl)
  for more functions like `n_observables(m)`, which returns
  the number of observables in `m`.