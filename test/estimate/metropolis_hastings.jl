using DSGE
using HDF5, Test
using DataFrames

path = dirname(@__FILE__)

# Set up model for testing
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")))
m = AnSchorfheide(custom_settings = custom_settings, testing = true)


# Read in the data, mode, and hessian
mode = load("$path/../reference/hessian.jld2","paramsmode")
data = load("$path/../reference/hessian.jld2","data")'

# Read in the covariance matrix for Metropolis-Hastings and reference parameter draws
hessian_inv, ref_draws, ref_cov =
    h5open("$path/../reference/metropolis_hastings.h5", "r") do file
        read(file, "hessian_inv"),
        read(file, "mhparams"),
        read(file, "ref_cov")
    end

# Set up and run metropolis-hastings
DSGE.update!(m, mode)
prop_cov = DegenerateMvNormal(mode, hessian_inv)

metropolis_hastings(prop_cov, m, data, .01, .09, verbose=:none)
compute_parameter_covariance(m)

# Read in the parameter draws and covariance just generated from estimate.
test_draws = h5open(rawpath(m, "estimate", "mhsave.h5"), "r") do file
    read(file, "mhparams")
end
test_cov = h5open(workpath(m, "estimate", "parameter_covariance.h5"), "r") do file
    read(file, "mhcov")
end

# Test that the parameter draws and covariance matrices are equal
@testset "Check equality of parameter draws and cov matrices in MH" begin
    @test @test_matrix_approx_eq ref_draws test_draws
    @test @test_matrix_approx_eq ref_cov test_cov
end

nothing
