"""
- `T, R, C`: transition equation matrices
- `Z, D`: observation equation matrices
- `Z_pseudo, D_pseudo`: matrices mapping states to psuedo-observables                                         
- `forecast_horizons`: select number of quarters ahead to forecast
- `variables_to_forecast`: select which forecasts to output ("States","Observables" and/or "Psuedo Observables") 
- `shocks_distribution`: specify either joint distribution from which to draw shocks
- `z`: state vector at time T
"""
function computeForecast(T::Array{Float64,2}, R::Array{Float64,2}, C::Array{Float64}, 
                         Z::Array{Float64,2}, D::Array{Float64},          
                         Z_pseudo::Array{Float64,2}, D_pseudo::Array{Float64},
                         forecast_horizons::Int64,
                         variables_to_forecast::Array,
                         shocks_distribution::Union{Distribution,Array{Float64,2}},
                         z::Array{Float64,1})
                                    
    nshocks=size(R,2)
    nstates=size(z,1)
    
    states=zeros(forecast_horizons,nstates)
    
    z_old=z
    
    for forecast_horizon=1:forecast_horizons

        # Step 1: Draw shocks
        if isa(shocks_distribution, Distribution)
            shocks = rand(shocks_distribution)
        elseif isa(shocks_distribution, Array{S,2})
            shocks = DistributionOfShocks
        end
        
        # Step 2: Forecast states
        z_new = C+T*z_old+R*shocks
        states[forecast_horizon,:]=z_new
        z_old=z_new
    end
    
    # Step 3: Calculate linear combinations of states (observables and pseudo-observables) as requested by user
    observables=Z*states'+repmat(D,1,forecast_horizons)
    psuedo_observables=Z_pseudo*states'+repmat(D_pseudo,1,forecast_horizons)
    
    # return a dictionary which houses all forecasts
    forecast=Dict("states"=>states, "observables"=>observables, "pseudo_observables"=>psuedo_observables)
    return forecast
    
end
