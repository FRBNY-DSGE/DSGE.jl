using HDF5
using Base.Test
using DSGE

include("../util.jl")
path = dirname(@__FILE__)

h5 = h5open("$path/solve.h5")
TTT_expected = read(h5, "TTT")
CCC_expected = reshape(read(h5, "CCC"), 78, 1)
RRR_expected = read(h5, "RRR")
close(h5)


model = Model990()
TTT, RRR, CCC = solve(model)
@test test_matrix_eq(TTT_expected, TTT)
@test test_matrix_eq(RRR_expected, RRR)
@test test_matrix_eq(CCC_expected, CCC)
