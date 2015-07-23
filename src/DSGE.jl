module DSGE

### HOW TO RUN THIS CODE
#
# From any directory:
#
#   using DSGE
#   using .M990
#   model = Model()


include("init/DistributionsExt.jl")
include("init/FinancialFrictionsFunctions.jl")
include("solve/Gensys.jl")

using Distributions, Roots.fzero, MATLAB
using .DistributionsExt, .FinancialFrictionsFunctions, .Gensys

import Base: convert, promote_rule, log, exp, start, next, done

export
    # core.jl
    AbstractModel, Param, toreal, tomodel, Parameters, prior, ModelInds, makedict,

    # solve/
    solve,

    # estimate/
    dlyap, kalcvf2NaN, kalsmth_k93, likelihood, posterior,

    # models/
    Parameters990, Model990, model_specifications, eqcond, measurement

include("core.jl")

include("solve/solve.jl")

include("estimate/kalman.jl")
include("estimate/posterior.jl")

include("models/m990/m990.jl")
include("models/m990/parameters.jl")
include("models/m990/modelinds.jl")

end
