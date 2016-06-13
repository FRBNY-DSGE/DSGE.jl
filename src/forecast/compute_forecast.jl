"""
```
compute_forecast(T, R, C, Z, D, Z_pseudo, D_pseudo, forecast_horizons,
variables_to_forecast, shocks, z)
```

## Inputs
- `T, R, C`: transition equation matrices
- `Z, D`: observation equation matrices
- `Z_pseudo, D_pseudo`: matrices mapping states to psuedo-observables                                         
- `forecast_horizons`: select number of quarters ahead to forecast
- `variables_to_forecast`: select which forecasts to output ("States","Observables" and/or
    "Psuedo Observables") 
- `shocks`: one-period joint distribution (type `Distribution`) from which to draw shocks or
    matrix of drawn shocks (size forecast_horizons x nshocks)
- `z`: state vector at time T

## Outputs
- dictionary of forecast outputs, with keys `:states`, `:observables`, and
    `:pseudo_observables`
"""
function compute_forecast(T::Array{Float64,2}, R::Array{Float64,2}, C::Array{Float64}, 
                         Z::Array{Float64,2}, D::Array{Float64},          
                         Z_pseudo::Array{Float64,2}, D_pseudo::Array{Float64},
                         forecast_horizons::Int,
                         variables_to_forecast::Array,
                         shocks::Array{Float64,2},
                         z::Array{Float64,1})

    if forecast_horizons <= 0
        throw(DomainError())
    end
                                    
    # Setup
    nshocks      = size(R,2)
    nstates      = size(T,2)
    nobservables = size(Z,1)
    npseudo      = size(Z_pseudo,1)
    states       = zeros(forecast_horizons,nstates)
    
    # Define our iteration function
    iterate(z_t1, ϵ_t) = C + T*z_t1 + R*ϵ_t

    # Iterate first period
    states[1, :] = iterate(z, shocks[1, :]')
    
    # Iterate remaining periods
    for t in 2:forecast_horizons
        states[t,:] = iterate(states[t-1,:]', shocks[t, :]')
    end

    # Apply observation and pseudo-observation equations
    observables        = D        .+ Z * states'
    pseudo_observables = D_pseudo .+ Z * states'
    
    # return a dictionary which houses all forecasts
    Dict{Symbol,Array{Float64}}(:states => states,
                                :observables => observables,
                                :pseudo_observables =>pseudo_observables)
end

# Utility method to actually draw shocks
function compute_forecast(T::Array{Float64,2}, R::Array{Float64,2}, C::Array{Float64}, 
                         Z::Array{Float64,2}, D::Array{Float64},          
                         Z_pseudo::Array{Float64,2}, D_pseudo::Array{Float64},
                         forecast_horizons::Int,
                         variables_to_forecast::Array,
                         dist::Distribution,
                         z::Array{Float64,1})

    if forecast_horizons <= 0
        throw(DomainError())
    end

    nshocks = size(R,2)
    shocks = zeros(forecast_horizons, nshocks)

    for t in 1:forecast_horizons
        shocks[t,:] = rand(dist)
    end

    compute_forecast(T, R, C, Z, D, Z_pseudo, D_pseudo, forecast_horizons,
        variables_to_forecast, shocks, z)
end
