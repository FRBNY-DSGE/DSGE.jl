using HDF5
using Base.Test
using DSGE

include("../util.jl")
path = dirname(@__FILE__)

h5 = h5open("$path/../reference/solve.h5")
TTT_expected = read(h5, "TTT")
CCC_expected = reshape(read(h5, "CCC"), 78, 1)
RRR_expected = read(h5, "RRR")
close(h5)


model = Model990()
TTT, RRR, CCC = solve(model)
@test_matrix_approx_eq TTT_expected TTT
@test_matrix_approx_eq RRR_expected RRR
@test_matrix_approx_eq CCC_expected CCC
