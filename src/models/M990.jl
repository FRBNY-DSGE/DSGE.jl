module M990

using DSGE.AbstractModel

export Parameters990, Model

include("m990/spec.jl")
include("m990/parameters.jl")
include("m990/modelinds.jl")
include("m990/eqcond.jl")

function Model()
    Θ = Parameters990()
    I = ModelInds(make_endo, make_exo, make_ex, make_eq)
    return AbstractModel.Model("990", Θ, I, eqcond)
end

end
