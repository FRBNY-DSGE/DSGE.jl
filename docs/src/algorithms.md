# Standard Algorithms

```@meta
CurrentModule = DSGE
```

## Solving the Model

```@docs
gensys
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

## State Space Filters
```@autodocs
Modules = [DSGE]
Pages   = ["kalman.jl"]
Order   = [:function, :type]
```

## State Space Smoothers
```@autodocs
Modules = [DSGE]
Pages   = ["smoothers.jl"]
Order   = [:function, :type]
```
