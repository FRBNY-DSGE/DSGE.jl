using Base: Test
using MATLAB

using DSGE
using DSGE: DistributionsExt
include("../util.jl")
path = dirname(@__FILE__)



mf = MatFile("$path/metropolis_hastings.mat")
mode = get_variable(mf, "params")
hessian = get_variable(mf, "hessian")
YY = get_variable(mf, "YYall")

σ = get_variable(mf, "sigscale")
Σ = get_variable(mf, "s")
rank = get_variable(mf, "sigpropdim")
logdet = get_variable(mf, "sigproplndet")

randvecs = get_variable(mf, "randvecs")
randvals = get_variable(mf, "randvals")
close(mf)

model = Model990()
cc0 = 0.01
cc = 0.09

if false
    propdist = DegenerateMvNormal(mode, σ, Σ, hessian, rank, logdet)
else
    propdist = proposal_distribution(mode, hessian)
end
para_sim = metropolis_hastings(propdist, model, YY, cc0, cc, randvecs, randvals)
mean_para_sim = mean(para_sim, 1)
