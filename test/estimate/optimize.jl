using DSGE
using HDF5
path = dirname(@__FILE__)
include("../util.jl")

# Test in model optimization
m = AnSchorfheide()
m.testing = true

# file = h5open("$path/../reference/optimize.h5","r")
# x0 = read(file, "params")
# data = read(file, "data")
# minimizer = read(file, "minimum")
# minimum = read(file, "f_minimum")
# H_expected = read(file, "H")
# close(file)

# # See src/estimate/estimate.jl
# update!(m, x0)
n_iterations = 3

# todo - remove
x0 = [p.value for p in m.parameters]

@time out, H = optimize!(m, data; iterations=n_iterations)

# write test results - TEMP
file = h5open("$path/../reference/optimize.h5","w")
file["params"] = x0
file["data"] = data
file["minimum"] = out.minimizer
file["f_minimum"] = out.minimum
file["H"] = H
close(file)

# @test_matrix_approx_eq minimizer out.minimizer
# @test_approx_eq_eps minimum out.minimum 5e-7
# @test_matrix_approx_eq H_expected H

nothing
