using DSGE
using HDF5, Base.Test
using DataFrames

include("../util.jl")
path = dirname(@__FILE__)

# Set up model for testing
m = Model990()
m.testing=true

# Read in the covariance matrix for Metropolis-Hastings and reference parameter draws
file = h5open("$path/../reference/metropolis_hastings.h5","r")
propdist_cov  = read(file, "propcov")
ref_draws     = read(file, "ref_draws")
ref_cov       = read(file, "ref_cov")
close(file)

# Read in the data
df = load_data(m; try_disk=true, verbose=:none)

# Set up and run metropolis-hastings
specify_mode!(m, inpath(m, "user", "paramsmode.h5"), verbose=:none)
specify_hessian(m, inpath(m, "user", "hessian.h5"), verbose=:none)

# Note Distributions (used in earlier test) also exports `estimate`.
DSGE.estimate(m, df; verbose=:none, proposal_covariance = propdist_cov, method = :MH)

# Read in the parameter draws and covariance just generated from estimate.
test_draws = h5open(rawpath(m, "estimate", "mhsave.h5"), "r") do file
    read(file, "mhparams")
end
test_cov = h5open(workpath(m, "estimate", "parameter_covariance.h5"), "r") do file
    read(file, "mhcov")
end

return

# Test that the parameter draws and covariance matrices are equal
@test_matrix_approx_eq ref_draws test_draws
@test_matrix_approx_eq ref_cov test_cov

# Test that the fixed parameters are all fixed
for fixed_param in [:δ, :λ_w, :ϵ_w, :ϵ_p, :g_star]
    @test test_draws[:,m.keys[fixed_param]] == fill(Float32(m[fixed_param].value), 100)
end

# Make sure that compute_moments runs appropriately
compute_moments(m, verbose=:none)

nothing
