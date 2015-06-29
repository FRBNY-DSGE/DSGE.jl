module M990

using DSGE.AbstractModel

export Parameters990, ModelInds990, Model990

include("m990/spec.jl")
include("m990/parameters.jl")
include("m990/modelinds.jl")
include("m990/eqcond.jl")

function Model990()
    Θ = Parameters990()
    I = ModelInds990()
    return AbstractModel.Model("990", Θ, I, eqcond)
end

end
