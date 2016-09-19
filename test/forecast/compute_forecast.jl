# include necessary packages
using DSGE
using Distributions

# Prepare system matrices, for the following hypothetical system
nstates      = 2
nequations   = 2
nshocks      = 1
nobservables = 1
npseudos     = 2
T        = eye(nequations,nstates)
R        = ones(nequations,nshocks)
C        = zeros(nequations)
Z        = ones(nobservables,nstates)
D        = zeros(nobservables)
Z_pseudo = ones(npseudos,nstates)
D_pseudo = zeros(npseudos)
Q        = eye(nshocks,nshocks)

# load state vector at time T
z_end = zeros(2)

# specify forecast horizons
forecast_horizons = 3

# define shock distribution
dist = DSGE.DegenerateMvNormal(zeros(nshocks), Q)

# test invocation supplying distribution
forecast_dist = DSGE.compute_forecast(T, R, C, Z, D, Z_pseudo, D_pseudo,
    forecast_horizons, dist, z_end)

# test invocation supplying shocks
shocks = rand(nshocks, forecast_horizons)
forecast_shocks = DSGE.compute_forecast(T, R, C, Z, D, Z_pseudo, D_pseudo,
    forecast_horizons, shocks, z_end)

nothing
