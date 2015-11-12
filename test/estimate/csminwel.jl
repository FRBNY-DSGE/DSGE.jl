using DSGE
using HDF5

include("../util.jl")
path = dirname(@__FILE__)

# Test in generic case
# Global minimum at (a, a^2)
a = 1
b = 100
function rosenbrock(x::Vector)
    return (a-x[1])^2.0 + b*(x[2]-x[1]^2.0)^2.0
end

# this is the actual gradient of the Rosenbrock function
function rosenbrock_grad(x::Array)
    dr = similar(x)
    dr[1] = -2*(a-x[1]) - 2b*(2x[1])*(x[2]-x[1]^2.0)
    dr[2] = 2b*(x[2]-x[1]^2.0)
    badg = false
    return dr, badg
end

# A really bad guess
x_init = [10.0, -9.0]
x_expected = [a, a^2.0]

res_real_grad, _ = csminwel(rosenbrock, rosenbrock_grad, [10.0, -9.0])
@test_approx_eq x_expected res_real_grad.minimum
res_numeric_grad, _ = csminwel(rosenbrock, [10.0, -9.0])
@test_approx_eq_eps x_expected res_numeric_grad.minimum 1e-8

# Test in model optimization
model = Model990()
model.testing=true

params = h5open(inpath(model, "user", "mode_in.h5")) do file
    read(file, "params")
end
YY = h5open(inpath(model, "data", "data_REF.h5")) do file
    read(file, "YY")
end
file = h5open("$path/../reference/csminwel_out.h5","r")
minimum_ = read(file, "minimum")
f_minimum = read(file, "f_minimum")
H_expected = read(file, "H_expected")
close(file)

# See src/estimate/estimate.jl
update!(model, params)
nit = 5
crit = 1e-10

out, H = optimize!(model, YY; ftol=crit, iterations=nit)

@test_matrix_approx_eq minimum_ out.minimum
@test_approx_eq f_minimum out.f_minimum
@test_matrix_approx_eq H_expected H

nothing
