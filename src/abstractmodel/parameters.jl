# The abstract Parameters type is the supertype of all model-specific ParametersXXX types.
# All concrete types have both Param (parameters) and Float64 (steady-state values) fields. See Parameters990 for an example.
abstract Parameters



# TODO: calculated logpdf values for inverse gamma don't correspond to priodens.m results
# Calculate (log of) joint density of Θ
function logprior(Θ::Parameters)
    sum = 0.0
    for i = 1:length(names(Θ))
        α = getfield(Θ, i)
        if !isa(α, Param)
            break
        else
            curr = logpdf(α.prior, α.value)
            #println(string(i, "   ", curr))
            sum += curr
        end
    end
    return sum
end
