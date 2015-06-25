module DistributionsExt
# This file defines additional functions to return objects of type Distribution. This is necessary because the original
# Matlab code specifies prior distributions wrt mean and SD (for beta and gamma-distributed parameters) and ν and σ (for
# inverse gamma-distributed parameters). Note that these functions are NOT new methods for the Distributions.Beta, etc.
# functions, but rather new functions with the same names.
# See also https://github.com/JuliaLang/julia/issues/9547

using Distributions

# Define PointMass distribution for fixed parameters.
type PointMass <: Distribution{Univariate, Continuous}
    μ::Float64
end

function Distributions.pdf(d::PointMass, x::Real)
    x == d.μ ? 1.0 : 0.0
end

# Given μ and σ, calculate α and β, return distribution
function Beta(μ::Real, σ::Real)
    α = (1-μ) * μ^2 / σ^2 - μ
    β = α * (1/μ - 1)
    return Distributions.Beta(α, β)
end

function Gamma(μ::Real, σ::Real)
    β = σ^2 / μ
    α = μ / β
    return Distributions.Gamma(α, β)
end

# ν and σ parameterize the scaled inverse chi-squared distribution, a different parameterization of the inverse gamma
# distribution.
function InverseGamma(σ::Real, ν::Real)
    α = ν / 2
    β = ν * σ^2 / 2
    return Distributions.InverseGamma(α, β)
end    

end
