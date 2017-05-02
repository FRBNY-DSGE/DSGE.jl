using DSGE
using HDF5

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

res_real_grad, _ = csminwel(rosenbrock, rosenbrock_grad, x_init)
@test_approx_eq x_expected res_real_grad.minimizer
res_numeric_grad, _ = csminwel(rosenbrock, x_init)
@test_approx_eq_eps x_expected res_numeric_grad.minimizer 1e-8

nothing
