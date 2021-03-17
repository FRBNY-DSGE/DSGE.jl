# [Solving the Model](@id solving-dsge-doc)

## The `gensys` routine

The DSGE model is written down in its canonical representation:

``` math
\Gamma_0 s_t = \Gamma_1 s_{t-1} + C + \Psi \epsilon_t + \Pi \eta_t
```

where ``\Gamma_0``, ``\Gamma_1``, ``C``, ``\Psi``, and ``\Pi`` are matrices of
coefficients for ``s_t`` (states at time ``t``), ``s_{t-1}`` (lagged states),
``\epsilon_t`` (exogenous shocks) and ``\eta_t`` (expectational shocks).

DSGE.jl solves the model to obtain its state-space form:

```math
\begin{aligned}
s_t &= T s_{t-1} + R \epsilon_t + C & \epsilon_t &\sim N(0, Q) & \mathrm{(transition)} \\
y_t &= Z s_t + D + u_t & u_t &\sim N(0, E) & \mathrm{(measurement)}
\end{aligned}
```

using the `gensys` routine of Chris Sims, introduced in
[this paper](http://sims.princeton.edu/yftp/gensys/LINRE3A.pdf).
This algorithm can be easily extended to (exogenous) regime switching.
For each regime ``i``, define the canonical matrices
``\Gamma_{0, i}``, ``\Gamma_{1, i}``, ``C_i``, ``\Psi_i``, and ``\Pi_i``.
Calling `gensys` on each regime and constructing the relevant easurement matrix
yields the matrices ``T_i``, ``R_i``,
``C_i``, ``Q_i``, ``Z_i``, ``D_i``, and ``E_i``, which define the state-space form
of a DSGE with multiple regimes.


We provide a
standalone native Julia implementation of the routine ([`gensys`](@ref)) as well
as a wrapper for `AbstractDSGEModel` subtypes ([`solve`](@ref)). When the Gensys.jl
package becomes ready for use, we intend to deprecate our `gensys` code and
substitute the `gensysdt` method for our code.

## [Regime-Switching](@id solveregswitch)
We allow solving the model with regime-switching in two cases.

1. Exogenous and unanticipated regimes
2. Alternative policies (permanent and temporary)

The first is straightforward to implement because
the regimes are exogenous and unanticipated. We simply need to specify
the equilibrium conditions in each regime, run `gensys` for each regime,
and return multiple transition equations. The required steps
to solve a model with exogenous and unanticipated regime-switching are

1. Write `eqcond` to accept a second argument specifying the regime, e.g. `eqcond(m::MyDSGEModel, reg::Int)`.
To allow no regime-switching, we recommend also writing the wrapper `eqcond(m) = eqcond(m, 1)`
or whatever default regime is desired.

2. Add additional keyword arguments to `measurement` so that it is defined as
```
function measurement(m::MyDSGEModel, TTT::AbstractMatrix{T}, RRR::AbstractMatrix{T}, CCC::AbstractVector{T};
                     reg::Int = 1, TTTs::Vector{<: AbstractMatrix{T}} = Matrix{T}[],
                     CCCs::Vector{<: AbstractVector{T}} = Vector{T}[],
                     information_set::UnitRange = reg:reg,
                     memo::Union{ForwardMultipleExpectationsMemo, Nothing} = nothing) where {T <: Real}
```
The type assertions for the arrays and keywords are not strictly necessary but are advised. Alternatively,
the user could simply set
```
function measurement(m::MyDSGEModel, TTT::AbstractMatrix{T}, RRR::AbstractMatrix{T}, CCC::AbstractVector{T};
                     kwargs...) where {T <: Real}
```
to avoid problems with forgetting certain kwargs.

3. Add the settings and indices required to ensure regime-switching is properly handled.
See [Regime-Switching Forecasts](@ref) for guidance on how to do this.

For 1 and 2, we recommend conferring with the implementation of regime-switching for Model 1002.
See the [equilibrium conditions here](https://github.com/FRBNY-DSGE/DSGE.jl/blob/main/src/models/representative/m1002/eqcond.jl)
and the [measurement equation here](https://github.com/FRBNY-DSGE/DSGE.jl/blob/main/src/models/representative/m1002/measurement.jl).

Temporary alternative policies are slightly more complicated. They leverage the same machinery
as exogenous and unanticipated regime-switching, so steps 1 and 2 above are required still.
But in addition, we need to specify in which regimes policies should be different
than those specified by `eqcond`, which assumes all policies are permanent ones.
Once we have specified these policies, we use the algorithm from
[Calgiarini and Kulish](https://www.mitpressjournals.org/doi/pdf/10.1162/REST_a_00240)
to calculate the rational expectations solution to a model with ``predicted structural changes'',
which allows us to specify temporary alternative policies like a temporary ZLB or
a temporary switch to average inflation targeting. The function implementing this algorithm
is `gensys2`.

```@docs
DSGE.gensys2
```

See [Alternative Policies](@ref)
and the regime-switching [example script](https://github.com/FRBNY-DSGE/DSGE.jl/blob/main/examples/regime_switching.jl)
for guidance on how to use temporary alternative policies.


```@docs
DSGE.solve
```

## The `System` Type
The `solve` function only returns the ``TTT``, ``RRR``, and ``CCC`` matrices.
To obtain the other matrices in the state-space representation of the DSGE,
we provide a wrapper function `compute_system`, which returns an object of
type `System`, whose fields are also special types we have implemented to
faciltiate use of the state-space form.

```@docs
System
Transition
Measurement
PseudoMeasurement
compute_system
```

## The `RegimeSwitchingSystem` Type
When there are multiple (exogenous) regimes in the state-space representation
of the DSGE, it is useful to treat the system as one type rather than
a vector of individual `System` objects. The `RegimeSwitchingSystem` type
implements this idea, and via multiple dispatch, most commands that work
on a `System` object also work on a `RegimeSwitchingSystem` object.

```@docs
RegimeSwitchingSystem
```

To get the number of regimes, one can call `n_regimes(regime_switch_system)` or
`regime_switch_system[:regimes]`.
Like a `System` object, the fields of a `RegimeSwitchingSystem` can be accessed by

```
regime_switch_system[:transitions]
regime_switch_system[:measurements]
regime_switch_system[:pseudo_measurements]
```

A `System` can be formed from a specific regime by calling

```
sys_reg2 = regime_switch_system[2]
```

and specific matrices/types in a given regime can be accessed by

```
regime_switch_system[2, :TTT]
regime_switch_system[2, :transition]
regime_switch_system[2, :measurement]
regime_switch_system[2, :pseudo_measurement]
```

## [The `TimeVaryingInformationSetSystem` Type](@id tvistype)
This type implements a state space system with [time-varying information sets](@ref tvis).
We alter two parts of the `RegimeSwitchingSystem`.

First, the `transitions` field
is now a `Vector{Vector{Transition}}`` because we may need to calculate the transition equations
under different information sets. For example, the time-varying transition equations implied by expecting
a temporary ZLB to 2023:Q4 in 2020:Q4 will be different than those implied by expecting a temporary ZLB to 2024:Q4.
Thus, if the Federal Reserve announced a ZLB to 2023:Q4 in 2020:Q4, and then in 2021:Q1, they announced
an extension of the ZLB to 2024:Q4, then the measurement equation needs to calculate
two sets of time-varying transition equations.

Second, we add two new fields, `information_set` and
`select`. The first specifies which regimes are part of the information set in a given regime.
The second specifies which transition equation to use in each regime. Continuing
the example of the change in the ZLB length, the regime in the `select` field for 2020:Q4
would point to the set of transition equations associated with a temporary ZLB to 20203:Q4,
while the regime in `select` for 2021:Q1 would point
to the set of transition equations associated with a temporary ZLB to 20204:Q4.

For more guidance, see [Time-Varying Information Sets](@ref tvis) and
the example [tvis_system.jl](https://github.com/FRBNY-DSGE/DSGE.jl), which
shows how to work with this type.
