using DSGE
using HDF5, Base.Test
using DataFrames

include("../util.jl")
path = dirname(@__FILE__)

# Set up model for testing
m = Schorf()
m.testing=true

# Read in the data
df = load_data(m; try_disk=true, verbose=:none)

# Set up and run metropolis-hastings
specify_mode!(m, inpath(m, "user", "paramsmode.h5"), verbose=:none)
specify_hessian(m, inpath(m, "user", "hessian.h5"), verbose=:none)

# Note Distributions (used in earlier test) also exports `estimate`.
DSGE.estimate(m, df; verbose=:none, proposal_covariance = propdist_cov)

# Testing of systematic resampling step.
file = h5open("$path/../reference/degen_dist.h5", "w")
wtsim = read(file, "wtsim")
npart = read(file, "npart")
close(file)

ESS = 1/sum(wtsim.^2)
if (ESS < npart/2)
    error("systematic resampling not working")
end

# Make sure that compute_moments runs appropriately
compute_moments(m, verbose=:none)

nothing


