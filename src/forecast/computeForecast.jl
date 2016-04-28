function computeForecast{S<:AbstractFloat}(
                                           T::Matrix{S}, R::Matrix{S}, C::Matrix{S},                  # transition equation matrices
                                           Z::Matrix{S}, D::Matrix{S},                                # observation equation matrices
                                           Z_pseudo::Matrix{S}, D_pseudo::Matrix{S},                  # matrices mapping states to psuedo-observables
                                           forecast_horizons::Int64,                               # select number of quarters ahead to forecast
                                           variables_to_forecast::Array{any,1},                       # select which forecasts to output ("States","Observables" and/or "Psuedo Observables")
                                           shocks_distribution::Union{Distribution,Matrix{S}})        # specify either joint distribution from which to draw shocks
                                           #z::vector{})                                               # state vector at time T
                                    
#return true                                    

# Step 1: Get shocks ready
#if isa(shocks_distribution, Distribution);
#    shocks = rand(distribution_of_shocks,[nshocks, forecast_horizons]);
#elseif isa(shocks_distribution, Matrix);
#    shocks = DistributionOfShocks;
#end

# Step 2: Forecast states
#z_new = CCC+TTT*z+RRR*Shocks(t,:)';

# Step 3: If wish to output observables, map forecasted states to observables
# add observable forecasts to dictionary


# Step 4: If wish to output mapped states, map forecast states to observables
# add mapped states forecasts to dictionary

# Step 5: output dictionary containing all results
# return a custom type called "forecast" containing different variables forecasts as well as shocks which were used

end
