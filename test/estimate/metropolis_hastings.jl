using DSGE
using HDF5, Base.Test

include("../util.jl")
path = dirname(@__FILE__)

# Set up model for testing
model = Model990()
model.testing=true

# Read in the covariance matrix for Metropolis-Hastings and reference parameter draws
file = h5open("$path/../reference/metropolis_hastings.h5","r")
propdist_cov  = read(file, "propdist_cov")
ref_draws     = read(file, "ref_draws")
ref_cov       = read(file, "ref_cov")
close(file)

println("Calling estimate")
# Set up and run metropolis-hastings
DSGE.estimate(model; verbose=:high, proposal_covariance = propdist_cov)

# Read in the parameter draws and covariance just generated from estimate()
test_draws = h5open(rawpath(model, "estimate", "mh_save.h5"), "r") do file
    read(file, "parasim")
end
test_cov = h5open(workpath(model, "estimate", "parameter_covariance.h5"), "r") do file
    read(file, "param_covariance")
end

# Test that the fixed parameters are all fixed
for fixed_param in [:δ, :λ_w, :ϵ_w, :ϵ_p, :g_star]
    @test test_draws[:,model.keys[fixed_param]] == fill(@compat(Float32(model[fixed_param].value)), 100)
end

# Test that the parameter draws are equal
@test_matrix_approx_eq ref_draws test_draws

# Test that the covariance matrices are equal
@test_matrix_approx_eq ref_cov test_cov

# Make sure that compute_moments runs appropriately
compute_moments(model, verbose=:none)

nothing
