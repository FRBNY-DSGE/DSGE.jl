using DSGE, ModelConstructors
using FileIO, HDF5, Test
using DataFrames

path = dirname(@__FILE__)
writing_output = false

# Set up model for testing
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")))
m = AnSchorfheide(custom_settings = custom_settings, testing = true)


# Read in the data, mode, and hessian
mode = load("$path/../reference/hessian.jld2", "paramsmode")
data = Matrix{Float64}(load("$path/../reference/hessian.jld2", "data")')

# Read in the covariance matrix for Metropolis-Hastings and reference parameter draws
# TODO: check that new MH agrees with old MH. (read in metropolis_hastings.h5)
hessian_inv, ref_draws, ref_cov =
    h5open("$path/../reference/metropolis_hastings_test.h5", "r") do file
        read(file, "hessian_inv"),
        read(file, "mhparams"),
        read(file, "ref_cov")
    end

# Set up and run metropolis-hastings
DSGE.update!(m, mode)
prop_cov = DegenerateMvNormal(mode, hessian_inv)

metropolis_hastings(prop_cov, m, data, .01, .09; verbose=:none)
compute_parameter_covariance(m)

# Read in the parameter draws and covariance just generated from estimate.
test_draws = h5open(rawpath(m, "estimate", "mhsave.h5"), "r") do file
    read(file, "mhparams")
end
test_cov = h5open(workpath(m, "estimate", "parameter_covariance.h5"), "r") do file
    read(file, "mhcov")
end

if writing_output
    h5open("$path/../reference/metropolis_hastings_test.h5", "w") do file
        write(file, "hessian_inv", hessian_inv),
        write(file, "mhparams",    test_draws),
        write(file, "ref_cov",     test_cov)
    end
end

# Test that the parameter draws and covariance matrices are equal
@testset "Check equality of parameter draws and cov matrices in MH (1 block)" begin
    @test @test_matrix_approx_eq ref_draws test_draws
    @test @test_matrix_approx_eq ref_cov test_cov
end

# Set up and run metropolis-hastings with three blocks!
hessian_inv, ref_draws, ref_cov =
    h5open("$path/../reference/metropolis_hastings_test_3_blocks.h5", "r") do file
        read(file, "hessian_inv"),
        read(file, "mhparams"),
        read(file, "ref_cov")
    end

DSGE.update!(m, mode)
prop_cov = DegenerateMvNormal(mode, hessian_inv)
m <= Setting(:n_mh_blocks, 3)

metropolis_hastings(prop_cov, m, data, .01, .09; verbose=:none)
compute_parameter_covariance(m)

# Read in the parameter draws and covariance just generated from estimate.
test_draws = h5open(rawpath(m, "estimate", "mhsave.h5"), "r") do file
    read(file, "mhparams")
end
test_cov = h5open(workpath(m, "estimate", "parameter_covariance.h5"), "r") do file
    read(file, "mhcov")
end

if writing_output
    h5open("$path/../reference/metropolis_hastings_test_3_blocks.h5", "w") do file
        write(file, "hessian_inv", hessian_inv),
        write(file, "mhparams",    test_draws),
        write(file, "ref_cov",     test_cov)
    end
end

# Test that the parameter draws and covariance matrices are equal
@testset "Check equality of parameter draws and cov matrices in MH (3 blocks)" begin
    @test @test_matrix_approx_eq ref_draws test_draws
    @test @test_matrix_approx_eq ref_cov test_cov
end

nothing
