# The given fields define the entire model structure.
# We can then concisely pass around a Model object to the remaining steps of the model
#   (solve, estimate, and forecast).
type Model990 <: AbstractModel
    spec_vars::Dict{String, Any}
    Θ::Parameters
    I::ModelInds
    eqcond::Function
    measurement::Function
end





function Model990()
    ### Model-specific specifications
    spec = Dict{String, Any}()
    # Number of anticipated policy shocks
    spec["nant"] = 6
    # Padding for nant
    spec["nantpad"] = 20
    # Number of periods back we should start incorporating zero bound expectations
    # ZLB expectations should begin in 2008 Q4
    spec["antlags"] = 24


    
    Θ = Parameters990(spec)
    I = ModelInds(spec)
    return Model990(spec, Θ, I, eqcond, measurement)
end

description(model::Model990) = "M990"
