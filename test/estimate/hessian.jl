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

hessian, _ = hessian!(model, mode, YY)
@test test_matrix_eq(hessian_expected[1:max_free_ind, 1:max_free_ind], 
                     hessian[1:max_free_ind, 1:max_free_ind]; 
                     ϵ=1.0)

model.testing = false

# Test `hessizero` in context of Rosenbrock function
function rosenbrock(x::Vector)
    a = 100
    b = 1
    return a*(x[2]-x[1]^2.0)^2.0 + b*(1-x[1])^2.0
end

function rosenbrock_hessian(x::Vector)
    H = zeros(eltype(x), 2, 2)
    H[1,1] = 2.0 - 400.0 * x[2] + 1200.0 * x[1]^2.0
    H[1,2] = -400.0 * x[1]
    H[2,1] = -400.0 * x[1]
    H[2,2] = 200.0
    return H
end

# At the min, ensure no negatives in diagonal
x0 = [1.0, 1.0]
hessian_expected = rosenbrock_hessian(x0)
hessian, _ = DSGE.hessizero(rosenbrock, x0; check_neg_diag=true)
@test test_matrix_eq(hessian_expected, hessian; ϵ=1e-2)

# Not at the min (indeed, we are at max), ensure throws error
x1 = [1.0, 1.0]
rosenbrock_neg(x) = -rosenbrock(x)
@test_throws Exception hessian, _ = DSGE.hessizero(rosenbrock_neg, x1; check_neg_diag=true)

# Not at the min, check matches closed form
x1 = Vector[[0.5, 1.5], [-1.0, -1.0]]
for x in x1
    hessian_expected = rosenbrock_hessian(x)
    hessian, _ = DSGE.hessizero(rosenbrock, x)
    @test test_matrix_eq(hessian_expected, hessian; ϵ=1e-2)
end
