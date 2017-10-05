# Solving the Model

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
\begin{align*}
s_t &= T s_{t-1} + R \epsilon_t + C & \epsilon_t &\sim N(0, Q) & \mathrm{(transition)} \\
y_t &= Z s_t + D + u_t & u_t &\sim N(0, E) & \mathrm{(measurement)}
\end{align*}
```

using the `gensys` routine of Chris Sims, introduced in
[this paper](http://sims.princeton.edu/yftp/gensys/LINRE3A.pdf). We provide a
standalone native Julia implementation of the routine ([`gensys`](@ref)) as well
as a wrapper for `AbstractModel` subtypes ([`solve`](@ref)). When the Gensys.jl
package becomes ready for use, we intend to deprecate our `gensys` code and
substitute the `gensysdt` method for our code.

```@docs
DSGE.solve
```
