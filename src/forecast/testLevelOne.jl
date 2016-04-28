
# include necessary packages
using Distributions
using DSGE

# load system matrices
using HDF5
T=h5read("fromMatlab.h5","T")
R=h5read("fromMatlab.h5","R")
C=h5read("fromMatlab.h5","C")
Z=h5read("fromMatlab.h5","Z")
D=h5read("fromMatlab.h5","D")
Z_pseudo=h5read("fromMatlab.h5","Z_pseudo")
D_pseudo=h5read("fromMatlab.h5","D_pseudo")
Q=h5read("fromMatlab.h5","Q")

# load state vector at time T
z_end=h5read("fromMatlab.h5","z")

# specify forecast horizons
forecast_horizons=60;

# specify which variables to forecast 
variables_to_forecast=Any["States","Observables","Pseudo-Observables"]

# define either shock distribution
shocks_distribution=DSGE.DegenerateMvNormal(zeros(size(Q,1)),Q)

function computeForecast(T::Array{Float64,2}, R::Array{Float64,2}, C::Array{Float64,1},           # transition equation matrices
                         Z::Array{Float64,2}, D::Array{Float64,1},                                # observation equation matrices
                         Z_pseudo::Array{Float64,2}, D_pseudo::Array{Float64,1},                  # matrices mapping states to psuedo-observables
                         forecast_horizons::Int64,                                                # select number of quarters ahead to forecast
                         variables_to_forecast::Array,                                            # select which forecasts to output ("States","Observables" and/or "Psuedo Observables")
                         shocks_distribution::Union{Distribution,Array{Float64,2}},               # specify either joint distribution from which to draw shocks
                         z::Array{Float64,1})                                                     # state vector at time T
                                    
                                    
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

## test computeForecast.jl
forecast=computeForecast(T, R, C, Z, D, Z_pseudo, D_pseudo, forecast_horizons, variables_to_forecast, shocks_distribution, z_end)


