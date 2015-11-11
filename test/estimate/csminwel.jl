using DSGE
using HDF5

include("../util.jl")
path = dirname(@__FILE__)

# Test in generic case
function rosenbrock(x::Vector)
    a = 100
    b = 1
    return a*(x[2]-x[1]^2.0)^2.0 + b*(1-x[1])^2.0
end

# this is the actual gradient of the Rosenbrock function
function rosenbrock_grad(x::Array)
    dr = similar(x)
    dr[1] = -2*(1 - x[1]) - 8*105*x[1]*(x[2]-x[1]^2.0)^3.0
    dr[2] = 4*105*(x[2]-x[1]^2.0)^3.0
    badg = false
    return dr, badg
end

# A really bad guess
x_init = [10.0, -9.0]
x_expected = [1.0, 1.0]

model = Model990()

res_real_grad, _ = csminwel(rosenbrock, rosenbrock_grad, [10.0, -9.0])
@test x_expected == res_real_grad.minimum
res_numeric_grad, _ = csminwel(rosenbrock, [10.0, -9.0], model=model)
@test x_expected == res_numeric_grad.minimum

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
x5 = read(file, "x5")
f5 = read(file, "f5")
f_calls = read(file, "f_calls")
g_calls = read(file, "g_calls")
H_expected = read(file, "H")
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

@test_matrix_approx_eq x5 out.minimum
@test_approx_eq f5 out.f_minimum
@test f_calls == out.f_calls
@test g_calls == out.g_calls
@test_matrix_approx_eq H_expected H

nothing
