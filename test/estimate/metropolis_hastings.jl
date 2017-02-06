using DSGE
using HDF5, Base.Test
using DataFrames

include("../util.jl")
path = dirname(@__FILE__)

# Set up model for testing
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")))
m = AnSchorfheide(custom_settings = custom_settings, testing = true)


# Read in the data, mode, and hessian
mode = h5read("$path/../reference/hessian.h5","paramsmode")
data = h5read("$path/../reference/hessian.h5","data")'

# Read in the covariance matrix for Metropolis-Hastings and reference parameter draws
hessian_inv   = h5read("$path/../reference/metropolis_hastings.h5", "hessian_inv")
ref_draws     = h5read("$path/../reference/metropolis_hastings.h5", "ref_draws")
ref_cov       = h5read("$path/../reference/metropolis_hastings.h5", "ref_cov")

# Set up and run metropolis-hastings
update!(m, mode)
prop_cov = DSGE.DegenerateMvNormal(mode, hessian_inv)

DSGE.metropolis_hastings(prop_cov, m, data, .01, .09, verbose=:none)
compute_parameter_covariance(m)

# Read in the parameter draws and covariance just generated from estimate.
test_draws = h5open(rawpath(m, "estimate", "mhsave.h5"), "r") do file
    read(file, "mhparams")
end
test_cov = h5open(workpath(m, "estimate", "parameter_covariance.h5"), "r") do file
    read(file, "mhcov")
end

# Test that the parameter draws and covariance matrices are equal
@test_matrix_approx_eq ref_draws test_draws
@test_matrix_approx_eq ref_cov test_cov

nothing