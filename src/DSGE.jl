### HOW TO RUN THIS CODE
#
# From any directory:
#
#   using DSGE
#   using .M990
#   model = Model()

module DSGE

include("init/DistributionsExt.jl")
include("init/FinancialFrictionsFunctions.jl")

using Distributions
using .DistributionsExt

import Base: convert, promote_rule, log, exp, start, next, done


export
    # functions
    solve, dlyap, likelihood,

    Parameters990, ModelInds, Model990,

    Param, Parameters, toreal, tomodel, logprior, model_specifications, eqcond, measurement, makedict,
    
    # estimate/kalman.jl
    kalcvf2NaN, kalsmth_k93


include("core.jl")

include("solve/Gensys.jl")
include("solve/solve.jl")

include("estimate/dlyap.jl")
include("estimate/kalman.jl")
include("estimate/likelihood.jl")



# m990
include("models/m990/m990.jl")
include("models/m990/parameters.jl")
include("models/m990/modelinds.jl")
include("models/m990/eqcond.jl")



end
