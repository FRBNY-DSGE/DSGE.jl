module AbstractModel

export
    # types
    Param, Parameters, ModelInds, Model

include("abstractmodel/param.jl")
include("abstractmodel/parameters.jl")
include("abstractmodel/modelinds.jl")
include("abstractmodel/model.jl")

end
