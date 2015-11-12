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

model = Model990()

res_real_grad, _ = csminwel(rosenbrock, rosenbrock_grad, [10.0, -9.0])
@test_approx_eq x_expected res_real_grad.minimum
res_numeric_grad, _ = csminwel(rosenbrock, [10.0, -9.0], model=model)
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

function posterior_min!{T<:AbstractFloat}(x::Vector{T})
    tomodel!(model,x)
    return -posterior(model, YY; catch_errors=true)[:post]
end

x0 = toreal(model.parameters)
H = 1e-4 * eye(DSGE.num_parameters(model))
nit = 5
crit = 1e-10

out, H = csminwel(posterior_min!, x0, H; model=model, ftol=crit, iterations=nit, show_trace=false, verbose=:none)
# h5 = h5open("/home/rcemjs04/.julia/v0.4/DSGE/test/reference/csminwel_out.h5","w")
# h5["minimum"] = out.minimum
# h5["f_minimum"] = out.f_minimum
# h5["H_expected"] = H
# close(h5)

@test_matrix_approx_eq minimum_ out.minimum
@test_approx_eq f_minimum out.f_minimum
@test_matrix_approx_eq H_expected H

nothing
