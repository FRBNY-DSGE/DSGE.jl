using Base: Test
using MATLAB

using DSGE
include("../util.jl")
path = dirname(@__FILE__)


mf = MatFile("$path/test_likelihood.mat")
YY = get_variable(mf, "YY")
lh_expected = get_variable(mf, "lnpy")
close(mf)


model = Model990()
lh = likelihood(model, YY)
