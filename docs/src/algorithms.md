# Standard Algorithms

```@meta
CurrentModule = DSGE
```

## Solving the Model

```@docs
gensys
gensys2
```

## [Optimization](@id algs-optimization)

```@docs
csminwel
optimize!
```

## Hessian Approximation

```@autodocs
Modules = [DSGE]
Pages   = ["hessian.jl", "hessizero.jl"]
Order   = [:function, :type]
```

## Sampling

```@docs
metropolis_hastings
```

## State Space Filters and Smoothers

See [StateSpaceRoutines.jl](https://github.com/FRBNY-DSGE/StateSpaceRoutines.jl).

## Sequential Monte Carlo

See [SMC.jl](https://github.com/FRBNY-DSGE/SMC.jl).
