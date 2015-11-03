using Base: Test
using DSGE, HDF5

include("../util.jl")
path = dirname(@__FILE__)

mode = h5open("$path/../reference/mode_in_optimized.h5") do file
    read(file, "mode")
end
YY = h5open("$path/../reference/YY.h5") do file
    read(file, "YY")
end
hessian_expected = h5open("$path/../reference/hessian_optimized.h5") do file
    read(file, "hessian")
end

model = Model990()

# Ad-hoc testing...
toggle_test_mode(model)
para_free      = [!θ.fixed for θ in model.parameters]
para_free_inds = find(para_free)
num_free_hessian_test = 4
para_free_inds = para_free_inds[1:num_free_hessian_test]
max_free_ind = max(para_free_inds...)

hessian, hessian_errors = hessian!(model, mode, YY; verbose=true)
@test test_matrix_eq(hessian_expected[1:max_free_ind, 1:max_free_ind], 
                     hessian[1:max_free_ind, 1:max_free_ind]; 
                     ϵ=1.0, verbose=true)

toggle_test_mode(model)
