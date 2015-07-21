### HOW TO RUN THIS CODE
#
# From any directory:
#
#   using DSGE
#   using .M990
#   model = Model()

module DSGE

export
    # functions
    solve, dlyap, likelihood,

    Parameters990, ModelInds, Model,

    Param, Parameters, toreal, tomodel, logprior, spec_vars, eqcond, measurement, makedict


include("abstractmodel/param.jl")
include("abstractmodel/parameters.jl")
include("abstractmodel/modelinds.jl")
include("abstractmodel/model.jl")

include("init/DistributionsExt.jl")
include("init/FinancialFrictionsFunctions.jl")
include("core.jl")

include("solve/Gensys.jl")
include("solve/solve.jl")

include("estimate/dlyap.jl")
include("estimate/Kalman.jl")
include("estimate/likelihood.jl")





include("models/m990/spec.jl")
include("models/m990/parameters.jl")
include("models/m990/modelinds.jl")
include("models/m990/eqcond.jl")
include("models/m990/measurement.jl")

function Model()
    spec = spec_vars()
    Θ = Parameters990(spec)
    I = ModelInds(spec)
    return Model("990", spec, Θ, I, eqcond, measurement)
end


end
