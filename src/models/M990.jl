module M990

using ..AbstractModel
import ..AbstractModel: Model

export Parameters990, ModelInds, Model, eqcond

include("m990/spec.jl")
include("m990/parameters.jl")
include("m990/modelinds.jl")
include("m990/eqcond.jl")

function AbstractModel.Model()
    Θ = Parameters990()
    I = ModelInds()
    return AbstractModel.Model("990", Θ, I, eqcond)
end

end
