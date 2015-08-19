module DSGE

using Compat
using Distributions, Roots.fzero, MATLAB

export
    # DSGE.jl
    savepath, inpath, outpath, tablepath, plotpath, logpath,

    # abstractdsgemodel.jl
    AbstractDSGEModel, Param, update!, toreal, tomodel, Parameters, tomodel!, prior, ModelInds, makedict,

    # solve/
    ordschur, gensys, solve,

    # estimate/
    dlyap!, kalcvf2NaN, kalsmth_k93, likelihood, posterior, posterior!, csminwel, hessizero!, estimate, proposal_distribution, metropolis_hastings,

    # models/
    steadystate!, Model990, model_specifications, eqcond, measurement

const savepath  = joinpath(dirname(@__FILE__), "../save/")
const inpath    = joinpath(savepath, "input_data/")
const outpath   = joinpath(savepath, "output_data/")
const tablepath = joinpath(savepath, "results/tables/")
const plotpath  = joinpath(savepath, "results/plots/")
const logpath   = joinpath(savepath, "logs/")

include("distributions_ext.jl")
include("abstractdsgemodel.jl")
include("parameters.jl")

if VERSION <= v"0.4-"
    include("solve/ordered_qz.jl")
end

include("solve/gensys.jl")
include("solve/solve.jl")

include("estimate/kalman.jl")
include("estimate/posterior.jl")
include("estimate/csminwel.jl")
include("estimate/hessian.jl")
include("estimate/estimate.jl")

include("models/m990/m990.jl")
include("models/m990/eqcond.jl")
include("models/m990/measurement.jl")

end
