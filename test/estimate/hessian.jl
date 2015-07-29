using Base: Test
using MATLAB

using DSGE
include("../util.jl")
path = dirname(@__FILE__)

mf = MatFile("$path/hessian.mat")
params = get_variable(mf, "params")
YY = get_variable(mf, "YY")
hessian_expected = get_variable(mf, "hessian")
close(mf)

model = Model990()

hessian, stoph = hessizero!(params, model, YY; noisy=true)
#@test test_matrix_eq(hessian_expected, hessian; noisy=true)
