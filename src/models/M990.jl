module M990

export Parameters990, ModelInds, Model

include("m990/spec.jl")
include("m990/parameters.jl")
include("m990/modelinds.jl")
include("m990/eqcond.jl")
include("m990/measurement.jl")

function Model()
    spec = spec_vars()
    Θ = Parameters990(spec)
    I = ModelInds(spec)
    return Model("990", spec, Θ, I, eqcond, measurement)
end

end
