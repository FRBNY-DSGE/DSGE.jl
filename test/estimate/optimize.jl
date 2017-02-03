using DSGE
using HDF5
path = dirname(@__FILE__)
include("../util.jl")

# Test in model optimization
m = AnSchorfheide()
m.testing = true

file = "$path/../reference/optimize.h5"
x0 = h5read(file, "params")
data = h5read(file, "data")'
minimizer = h5read(file, "minimizer")
minimum = h5read(file, "minimum")
H_expected = h5read(file, "H")

# See src/estimate/estimate.jl
update!(m, x0)
n_iterations = 3

x0 = Float64[p.value for p in m.parameters]

@time out, H = optimize!(m, data; iterations=n_iterations)

@test_matrix_approx_eq minimizer out.minimizer
@test_approx_eq_eps minimum out.minimum 5e-7
@test_matrix_approx_eq H_expected H

nothing
