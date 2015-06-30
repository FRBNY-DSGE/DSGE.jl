module AbstractModel

export
    Param, toreal, tomodel, # param.jl
    Parameters, logprior, # parameters.jl
    ModelInds, makedict, # modelinds.jl
    Model # model.jl

include("abstractmodel/param.jl")
include("abstractmodel/parameters.jl")
include("abstractmodel/modelinds.jl")
include("abstractmodel/model.jl")

end
