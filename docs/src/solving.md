# Solving the Model

## The `gensys` routine

The DSGE model is written down in its canonical representation:

```
Γ0 s_t = Γ1 s_{t-1} + C + Ψ ε_t + Π η_t
```

where Γ0, Γ1, C, Ψ, and Π are matrices of coefficients for s_t (states
at time t), s_{t-1} (lagged states), ε_t (exogenous shocks) and η_t
(expectational shocks).

DSGE.jl solves the model to obtain its state-space form using the
`gensys` routine of Chris Sims, introduced in [this
paper](http://sims.princeton.edu/yftp/gensys/LINRE3A.pdf). We provide
a standalone native Julia implementation of the routine ([`gensys`](@ref)) as
well as a wrapper for `AbstractModel` subtypes ([`solve`](@ref)):

When the Gensys.jl package becomes ready for use, we intend to
deprecate our `gensys` code and substitute the `gensysdt` method for
our code.

```@docs
DSGE.solve
```

