using MATLAB
using DSGE
path = dirname(@__FILE__)

mf = MatFile("$path/solve.mat")
TTT_expected = get_variable(mf, "TTT")
CCC_expected = reshape(get_variable(mf, "CCC"), 78, 1)
RRR_expected = get_variable(mf, "RRR")
close(mf)

model = Model990()
TTT, CCC, RRR = solve(model)
@test test_matrix_eq(TTT_expected, TTT)
@test test_matrix_eq(CCC_expected, CCC)
@test test_matrix_eq(RRR_expected, RRR)

