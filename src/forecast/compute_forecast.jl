"""
```
compute_forecast(T, R, C, Z, D, Z_pseudo, D_pseudo, forecast_horizons,
    shocks, z)
```

### Inputs

- `T`, `R`, `C`: transition equation matrices
- `Z`, `D`: observation equation matrices
- `Z_pseudo`, `D_pseudo`: matrices mapping states to pseudo-observables
- `forecast_horizons`: number of quarters ahead to forecast output
- `shocks`: joint distribution (type `Distribution`) from which to draw
  time-invariant shocks or matrix of drawn shocks (size `nshocks` x
  `forecast_horizons`)
- `z`: state vector at time `T`, i.e. at the beginning of the forecast

### Outputs

`compute_forecast` returns a dictionary of forecast outputs, with keys:

-`:states`
-`:observables`
-`:pseudo_observables`
-`:shocks`
"""
function compute_forecast{S<:AbstractFloat}(T::Matrix{S}, R::Matrix{S}, C::Vector{S}, 
                                            Z::Matrix{S}, D::Vector{S},
                                            Z_pseudo::Matrix{S}, D_pseudo::Vector{S},
                                            forecast_horizons::Int,
                                            shocks::Matrix{S},
                                            z::Vector{S})

    if forecast_horizons <= 0
        throw(DomainError())
    end
                                    
    # Setup
    nshocks      = size(R, 2)
    nstates      = size(T, 2)
    nobservables = size(Z, 1)
    npseudo      = size(Z_pseudo, 1)
    states       = zeros(nstates, forecast_horizons)
    
    # Define our iteration function
    iterate(z_t1, ϵ_t) = C + T*z_t1 + R*ϵ_t

    # Iterate first period
    states[:, 1] = iterate(z, shocks[:, 1])
    
    # Iterate remaining periods
    for t in 2:forecast_horizons
        states[:, t] = iterate(states[:, t-1], shocks[:, t])
    end

    # Apply observation and pseudo-observation equations
    observables        = D        .+ Z        * states
    pseudo_observables = D_pseudo .+ Z_pseudo * states
    
    # Return a dictionary of forecasts
    Dict{Symbol, Matrix{S}}(
        :states             => states,
        :observables        => observables,
        :pseudo_observables => pseudo_observables,
        :shocks             => shocks)
end

# Utility method to actually draw shocks
function compute_forecast{S<:AbstractFloat}(T::Matrix{S}, R::Matrix{S}, C::Vector{S}, 
                                            Z::Matrix{S}, D::Vector{S},          
                                            Z_pseudo::Matrix{S}, D_pseudo::Vector{S},
                                            forecast_horizons::Int,
                                            dist::Distribution,
                                            z::Vector{S})

    if forecast_horizons <= 0
        throw(DomainError())
    end

    nshocks = size(R, 2)
    shocks = zeros(nshocks, forecast_horizons)

    for t in 1:forecast_horizons
        shocks[:, t] = rand(dist)
    end

    compute_forecast(T, R, C, Z, D, Z_pseudo, D_pseudo, forecast_horizons,
        shocks, z)
end
