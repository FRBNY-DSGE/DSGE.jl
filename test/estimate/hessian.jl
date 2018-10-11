using DSGE
using HDF5, Test

path = dirname(@__FILE__)

# # Test hessian! in context of model
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")))
m = AnSchorfheide(custom_settings = custom_settings, testing = true)

# Setup paths

data = load("$path/../reference/hessian.jld2","data")'
mode = load("$path/../reference/hessian.jld2","paramsmode")
hessian_expected = load("$path/../reference/hessian.jld2","hessian")

# Test subset of hessian elements.
para_free      = [!θ.fixed for θ in m.parameters]
para_free_inds = findall(para_free)

max_free_ind = DSGE.n_hessian_test_params(m)
if max_free_ind < maximum(para_free_inds)
    para_free_inds = para_free_inds[1:max_free_ind]
end

hessian, _ = hessian!(m, mode, data; verbose=:none)

# The values here are liable to differ between machines, based on different
# architectures/multithreading. Hessian values tend to be large, so these larger tolerances
# are needed.
expect = hessian_expected[1:max_free_ind, 1:max_free_ind]
actual = hessian[1:max_free_ind, 1:max_free_ind]

@testset "Check Hessian calculation" begin
    @test @test_matrix_approx_eq_eps expect actual 0.1 3.0
end

m.testing = false

nothing
