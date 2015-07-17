module M990

using ..AbstractModel
import ..AbstractModel: Model

export Parameters990, ModelInds, Model

include("m990/spec.jl")
include("m990/parameters.jl")
include("m990/modelinds.jl")
include("m990/eqcond.jl")
include("m990/measurement.jl")

function AbstractModel.Model()
    spec = spec_vars()
    Θ = Parameters990(spec)
    I = ModelInds(spec)
    return AbstractModel.Model("990", spec, Θ, I, eqcond, measurement)
end

end
