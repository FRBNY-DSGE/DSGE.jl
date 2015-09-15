#using MATLAB
using HDF5

path = dirname(@__FILE__)

## mf = MatFile("$path/solve.mat")
## TTT_expected = get_variable(mf, "TTT")
## CCC_expected = reshape(get_variable(mf, "CCC"), 78, 1)
## RRR_expected = get_variable(mf, "RRR")
## close(mf)

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
