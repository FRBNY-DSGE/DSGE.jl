using Base: Test
using DSGE, HDF5

include("../util.jl")
path = dirname(@__FILE__)

h5 = h5open("$path/hessian.h5")
params = read(h5, "params")
YY = read(h5, "YY")
hessian_expected = read(h5, "hessian")
close(h5)

model = Model990()

@time hessian, stoph = hessizero!(model, params, YY; verbose=true)
@test test_matrix_eq(hessian_expected, hessian; ϵ=1.0, noisy=true)

