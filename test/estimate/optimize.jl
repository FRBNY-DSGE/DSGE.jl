using DSGE
using HDF5
path = dirname(@__FILE__)

# Test in model optimization
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")))
m = AnSchorfheide(custom_settings = custom_settings, testing = true)


file = "$path/../reference/optimize.jld2"
x0 = load(file, "params")
data = load(file, "data")'
minimizer = load(file, "minimizer")
minimum = load(file, "minimum")
H_expected = load(file, "H")

# See src/estimate/estimate.jl
DSGE.update!(m, x0)
n_iterations = 3

x0 = Float64[p.value for p in m.parameters]

out, H = optimize!(m, data; iterations=n_iterations)

# Re-generate test file
#=h5open("$path/../reference/optimize.h5", "w") do file
    file["params"] = params_test
    file["data"] = data_test
    file["minimizer"] = out.minimizer
    file["minimum"] = out.minimum
    file["H"] = H
end=#

@testset "Check optimize minimizers are the same" begin
    @test @test_matrix_approx_eq minimizer out.minimizer
    @test minimum ≈ out.minimum atol=5e-7
    @test @test_matrix_approx_eq H_expected H
end

nothing
