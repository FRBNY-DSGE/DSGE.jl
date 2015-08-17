module DSGE

using Compat
using Distributions, Roots.fzero, MATLAB

import Base: convert, promote_rule, log, exp, start, next, done, length

export
    # DSGE.jl
    savepath,

    # core.jl
    AbstractModel, Param, update!, toreal, tomodel, Parameters, prior, ModelInds, makedict,

    # solve/
    ordschur, gensys, solve,

    # estimate/
    dlyap!, kalcvf2NaN, kalsmth_k93, likelihood, posterior, posterior!, hessizero!, estimate,

    # models/
    steadystate!, Model990, model_specifications, eqcond, measurement

const savepath = joinpath(dirname(@__FILE__), "save/")

include("distributions_ext.jl")
include("core.jl")

include("solve/ordered_qz.jl")
include("solve/gensys.jl")
include("solve/solve.jl")

include("estimate/kalman.jl")
include("estimate/posterior.jl")
include("estimate/hessian.jl")
include("estimate/estimate.jl")

include("models/m990/m990.jl")
include("models/m990/eqcond.jl")
include("models/m990/measurement.jl")

# include("models/m990/parameters.jl")
# include("models/m990/modelinds.jl")

end
