using Base: Test
using DSGE, HDF5

include("../util.jl")
path = dirname(@__FILE__)

mode = h5open("$path/../reference/mode_in.h5") do file
    read(file, "params")
end
YY = h5open("$path/../reference/YY.h5") do file
    read(file, "YY")
end
hessian_expected = h5open("$path/../reference/hessian.h5") do file
    read(file, "hessian")
end

model = Model990()

@time hessian, stoph = hessizero!(model, mode, YY; verbose=true)
@test test_matrix_eq(hessian_expected, hessian; ϵ=1.0, verbose=true)

