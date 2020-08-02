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