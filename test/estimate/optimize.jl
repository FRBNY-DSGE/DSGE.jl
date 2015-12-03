using DSGE
using HDF5
path = dirname(@__FILE__)
include("../util.jl")

# Test in model optimization
m = Model990()
m.testing=true

file = h5open("$path/../reference/optimize.h5","r")
x0 = read(file, "params")
data = read(file, "data")
minimum_ = read(file, "minimum")
f_minimum = read(file, "f_minimum")
H_expected = read(file, "H")
close(file)

# See src/estimate/estimate.jl
update!(m, x0)
n_iterations = 3

@time out, H = optimize!(m, data; iterations=n_iterations)

@test_matrix_approx_eq minimum_ out.minimum
@test_approx_eq_eps f_minimum out.f_minimum 10e-8
@test_matrix_approx_eq H_expected H

nothing
