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
z0 = zeros(2)

# specify forecast horizons
horizons = 3

# define shock distribution
dist = DSGE.DegenerateMvNormal(zeros(nshocks), Q)

# test invocation supplying distribution
states, obs, pseudo, shocks = compute_forecast(T, R, C, Z, D, horizons, dist,
    z0, Z_pseudo, D_pseudo)

# test invocation supplying shocks
shocks = rand(nshocks, horizons)
states, obs, pseudo, shocks = compute_forecast(T, R, C, Z, D, horizons, shocks,
    z0, Z_pseudo, D_pseudo)

nothing
