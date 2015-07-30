module DSGE

include("init/DistributionsExt.jl")
include("init/FinancialFrictionsFunctions.jl")
include("solve/Gensys.jl")

using Distributions, Roots.fzero, MATLAB
using .DistributionsExt, .FinancialFrictionsFunctions, .Gensys

import Base: convert, promote_rule, log, exp, start, next, done, length

export
    # DSGE.jl
    savepath,

    # core.jl
    AbstractModel, Param, update!, toreal, tomodel, Parameters, prior, ModelInds, makedict,

    # solve/
    solve,

    # estimate/
    dlyap!, kalcvf2NaN, kalsmth_k93, likelihood, posterior, posterior!, hessizero!, estimate,

    # models/
    Parameters990, steadystate!, Model990, model_specifications, eqcond, measurement

const savepath = joinpath(dirname(@__FILE__), "save/")

include("core.jl")

include("solve/solve.jl")

include("estimate/kalman.jl")
include("estimate/posterior.jl")
include("estimate/hessian.jl")
include("estimate/estimate.jl")

include("models/m990/m990.jl")
include("models/m990/parameters.jl")
include("models/m990/modelinds.jl")

end
