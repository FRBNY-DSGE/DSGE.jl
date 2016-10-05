# DSGE.jl Documentation

## Solving the model
```@docs
solve(m::AbstractModel)
```

## Estimating the model
```@docs
DSGE.estimate(m::AbstractModel, df::DataFrame; verbose::Symbol=:low, proposal_covariance=Matrix())
DSGE.estimate(m::AbstractModel;verbose::Symbol=:low, proposal_covariance::Matrix=Matrix())
DSGE.estimate(m::AbstractModel, data::Matrix{Float64}; verbose::Symbol=:low, proposal_covariance::Matrix=Matrix())
metropolis_hastings{T<:AbstractFloat}(propdist::Distribution, m::AbstractModel, data::Matrix{T}, cc0::T, cc::T; verbose::Symbol=:low)
```


