using Base: Test
using MATLAB

using DSGE
using DSGE: M990
include("../util.jl")


mf = MatFile("test_likelihood.mat")
YY = get_variable(mf, "YY")
lh_expected = get_variable(mf, "lnpy")
close(mf)


model = Model()
lh = likelihood(model, YY)
