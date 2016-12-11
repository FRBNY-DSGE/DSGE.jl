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
minimizer = read(file, "minimum")
minimum = read(file, "f_minimum")
H_expected = read(file, "H")
close(file)

# See src/estimate/estimate.jl
update!(m, x0)
n_iterations = 3

@time out, H = optimize!(m, data; iterations=n_iterations)

@test_matrix_approx_eq minimizer out.minimizer
@test_approx_eq_eps minimum out.minimum 5e-7
@test_matrix_approx_eq H_expected H

nothing
