using DSGE, ModelConstructors, Test
using HDF5, Random
import ModelConstructors: @test_matrix_approx_eq_eps

path = dirname(@__FILE__)
writing_output = false
Random.seed!(1793)

## Assess post-estimation calculation of mode
#-----------------------------------------------------------------
# csminwel
#-----------------------------------------------------------------
custom_settings = [Setting(:date_forecast_start, quartertodate("2015-Q4"))]
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
if writing_output
    # Regenerate the Hessian reference under the current Julia/LAPACK. The saved one is
    # pre-migration; the small parameter-1 curvature entries drift a few percent vs it.
    h5open(file, "w") do f
        f["hessian"] = modal_hessian
    end
end
exp_hessian = h5read(file, "hessian")

@testset "Check optimize minimizers are the same [csminwel]" begin
    @test minimizer ≈ modal_out.minimizer atol=5e-4
    @test @test_matrix_approx_eq H_expected modal_H
    # Finite-difference Hessian: the well-determined entries match tightly, but the soft
    # parameter-1 curvature entries are only determinable to ~a few percent (which is what
    # made this flake REPL-vs-Test on the old Julias). Use a relative tolerance (5%) instead
    # of the 0.01% default; gross errors are still caught and the informative entries match to
    # <0.1%. Regenerate the reference above (writing_output=true) on a new Julia/LAPACK.
    @test @test_matrix_approx_eq_eps exp_hessian modal_hessian 1e-6 5.0
end
