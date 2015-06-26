include("parameters.jl")
include("modelinds.jl")

# The given fields define the entire model structure.
# We can then concisely pass around a Model object to the remaining steps of the model (solve, estimate, and forecast).
type Model
    spec::String
    Θ::Parameters
    I::ModelInds
    eqcond::Function
end

function Model(spec::String)
    include("../models/m$(spec)/spec.jl")
    include("../models/m$(spec)/parameters.jl")
    include("../models/m$(spec)/modelinds.jl")
    include("../models/m$(spec)/eqcond.jl")

    eval(parse("Θ = Parameters$(spec)()"))
    I = ModelInds(make_endo, make_exo, make_ex, make_eq)
    return Model(spec, Θ, I, eqcond)
end

