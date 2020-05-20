using DSGE, ModelConstructors, Test
using HDF5, Random
path = dirname(@__FILE__)
writing_output = false
Random.seed!(1793)

## Assess post-estimation calculation of mode
#-----------------------------------------------------------------
# csminwel
#-----------------------------------------------------------------
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")))
m = AnSchorfheide(custom_settings = custom_settings, testing = true)

# Load data
file = "$path/../reference/optimize.h5"
x0   = h5read(file, "params")
data = h5read(file, "data")'

# For regenerating test file
params_test = deepcopy(x0)
data_test   = Matrix{Float64}(data)

minimizer  = h5read(file, "minimizer")
minimum    = h5read(file, "minimum")
H_expected = h5read(file, "H")

# See src/estimate/estimate.jl
DSGE.update!(m, x0)
n_iterations = 3
m <= Setting(:optimization_iterations, n_iterations)
m <= Setting(:optimization_ftol, 1e-14)
m <= Setting(:optimization_xtol, 1e-32)
m <= Setting(:optimization_gtol, 1e-8)
m <= Setting(:optimization_attempts, 0)
m <= Setting(:optimization_step_size, .01)

modal_minimizer, modal_out, modal_H, modal_hessian =
    DSGE.calculate_mode(m, data_test, vec(params_test), :csminwel, mle = false,
                        verbose = :none, get_all_results = true,
                        save_results = false)

file = "$path/../reference/estimate_util.h5"
exp_hessian = h5read(file, "hessian")

@testset "Check optimize minimizers are the same [csminwel]" begin
    @test minimizer ≈ modal_out.minimizer atol=5e-4
    @test @test_matrix_approx_eq H_expected modal_H
    @test_broken @test_matrix_approx_eq exp_hessian modal_hessian # this works when ran in REPL but gives a different result in Test mode
end
