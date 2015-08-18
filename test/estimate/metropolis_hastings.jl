using Base: Test
using MATLAB
using HDF5

using DSGE
using DSGE: DistributionsExt
include("../util.jl")
path = dirname(@__FILE__)



mf = MatFile("$path/metropolis_hastings.mat")
mode = get_variable(mf, "params")
hessian = get_variable(mf, "hessian")
YY = get_variable(mf, "YYall")

randvecs = get_variable(mf, "randvecs")
randvals = get_variable(mf, "randvals")
close(mf)

mf2 = MatFile("$path/sigscale.mat")
σ = get_variable(mf2, "sigscale")
close(mf2)

model = Model990()
cc0 = 0.01
cc = 0.09

propdist = DegenerateMvNormal(mode, σ)

metropolis_hastings(propdist, model, YY, cc0, cc, randvecs, randvals)

