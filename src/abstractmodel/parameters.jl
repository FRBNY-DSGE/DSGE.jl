import Base: start, next, done
using Distributions: logpdf

# The abstract Parameters type is the supertype of all model-specific ParametersXXX types.
# All concrete types have both Param (parameters) and Float64 (steady-state values) fields.
# See Parameters990 for an example.
abstract Parameters



# Implement the iterator protocol for the Parameters type
# This will iterate over all Param fields (not steady-state values)
Base.start(Θ::Parameters) = 1
Base.next(Θ::Parameters, state::Int) = getfield(Θ, state), state+1
Base.done(Θ::Parameters, state::Int) = !isa(getfield(Θ, state), Param)



# TODO: calculated logpdf values for inverse gamma don't correspond to priodens.m results
# Calculate (log of) joint density of Θ
function logprior(Θ::Parameters)
    sum = 0.0
    for α = Θ
        curr = logpdf(α.prior, α.value)
        sum += curr
    end
    return sum
end
