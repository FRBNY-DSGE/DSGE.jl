
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

## test computeForecast.jl
forecast=computeForecast(T, R, C, Z, D, Z_pseudo, D_pseudo, forecast_horizons, variables_to_forecast, shocks_distribution, z_end)


