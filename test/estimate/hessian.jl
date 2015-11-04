using Base: Test
using DSGE, HDF5

include("../util.jl")
path = dirname(@__FILE__)

mode = h5open("$path/../reference/mode_in_optimized.h5") do file
    read(file, "mode")
end
YY = h5open("$path/../reference/data_REF.h5") do file
    read(file, "YY")
end
hessian_expected = h5open("$path/../reference/hessian_optimized.h5") do file
    read(file, "hessian")
end

# Test hessian! in context of model
model = Model990()

# Test subset of hessian elements.
model.testing = true

para_free      = [!θ.fixed for θ in model.parameters]
para_free_inds = find(para_free)

max_free_ind = DSGE.max_hessian_free_params(model)
if max_free_ind < maximum(para_free_inds)
    para_free_inds = para_free_inds[1:max_free_ind]
end

@time hessian, _ = hessian!(model, mode, YY; verbose=true)

expect = hessian_expected[1:max_free_ind, 1:max_free_ind]
actual = hessian[1:max_free_ind, 1:max_free_ind]
@test test_matrix_eq(expect, actual; ϵ_abs=1.0, ϵ_rel=1e-1)
@test_matrix_approx_eq_eps expect actual 1e-3 1e-1
# TODO calibrate above for when use_parallel_workers ≡ true

model.testing = false

nothing
