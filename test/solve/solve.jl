using MATLAB
using DSGE: M990

mf = MatFile("solve.mat")
TTT_expected = get_variable(mf, "TTT")
CCC_expected = reshape(get_variable(mf, "CCC"), 78, 1)
RRR_expected = get_variable(mf, "RRR")
close(mf)

model = Model()
TTT, CCC, RRR = solve(model)
@test test_matrix_eq(TTT_expected, TTT)
@test test_matrix_eq(CCC_expected, CCC)
@test test_matrix_eq(RRR_expected, RRR)

