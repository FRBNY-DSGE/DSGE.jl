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
