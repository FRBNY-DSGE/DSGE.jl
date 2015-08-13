# This file defines additional functions to return objects of type Distribution. This is
# necessary because the original Matlab code specifies prior distributions wrt mean and SD
# (for beta and gamma-distributed parameters) and ν and σ (for inverse gamma-distributed
# parameters).

# Define PointMass distribution for fixed parameters.
type PointMass <: Distribution{Univariate, Continuous}
    μ::Float64
end

Distributions.pdf(d::PointMass, x::Real) = (x == d.μ ? 1.0 : 0.0)
Distributions.mean(d::PointMass) = d.μ
Distributions.std(d::PointMass) = 0

# Given μ and σ, calculate α and β, return distribution
function BetaAlt(μ::Real, σ::Real)
    α = (1-μ) * μ^2 / σ^2 - μ
    β = α * (1/μ - 1)
    return Distributions.Beta(α, β)
end

function GammaAlt(μ::Real, σ::Real)
    β = σ^2 / μ
    α = μ / β
    return Distributions.Gamma(α, β)
end

# If x ~ RootInverseGamma(ν, τ²), then
#   x² ~ ScaledInverseChiSquared(ν, τ²)
#   x² ~ InverseGamma(ν/2, ντ²/2)

type RootInverseGamma <: Distribution{Univariate, Continuous}
    ν::Float64
    τ::Float64
end

Distributions.params(d::RootInverseGamma) = (d.ν, d.τ)

function Distributions.pdf(d::RootInverseGamma, x::Real)
    (ν, τ) = params(d)
    return 2 * (ν*τ^2/2)^(ν/2) * exp((-ν*τ^2)/(2x^2)) / gamma(ν/2) / x^(ν+1)
end

function Distributions.logpdf(d::RootInverseGamma, x::Real)
    (ν, τ) = params(d)
    return log(2) - log(gamma(ν/2)) + (ν/2)*log(ν*τ^2/2) - ((ν+1)/2)*log(x^2) - ν*τ^2/(2x^2)
end
