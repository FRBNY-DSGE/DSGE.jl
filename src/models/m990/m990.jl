# The given fields define the entire model structure.
# We can then concisely pass around a Model object to the remaining steps of the model
#   (solve, estimate, and forecast).
type Model990 <: AbstractModel
    spec::String
    spec_vars::Dict{String, Any}
    Θ::Parameters
    I::ModelInds
    eqcond::Function
    measurement::Function
end





function Model990()
    spec = spec_vars()
    Θ = Parameters990(spec)
    I = ModelInds(spec)
    return Model990("990", spec, Θ, I, eqcond, measurement)
end
