# This function initializes the model to a given specification and returns a Model object.
# It will eventually be called by the top-level code/main.jl file, but right now can be called with no other previous initialization:
#
#   m990 = init("990")

function init(spec::String)
    include("../abstractmodel/model.jl")

    include("../../models/m$(spec)/spec.jl")
    include("../../models/m$(spec)/parameters.jl")
    include("../../models/m$(spec)/modelinds.jl")
    include("../../models/m$(spec)/eqcond.jl")

    Θ = Parameters990()
    I = ModelInds(make_endo, make_exo, make_ex, make_eq)

    return Model(spec, Θ, I, eqcond)
end
