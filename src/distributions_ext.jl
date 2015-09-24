# This file defines additional functions to return objects of type Distribution. This is
# necessary because the original Matlab code specifies prior distributions wrt mean and SD
# (for beta and gamma-distributed parameters) and ν and σ (for inverse gamma-distributed
# parameters). Note these functions are NOT new methods for the Distributions.Beta, etc.
# functions, but rather new functions with the same names.

using Distributions, Compat
import Distributions: params, mean, std, pdf, logpdf, rand
import Base: length

export PointMass, Beta, Gamma, RootInverseGamma, DegenerateMvNormal
    
# Define PointMass distribution for fixed parameters.
type PointMass <: Distribution{Univariate, Continuous}
    μ::Float64
end

Distributions.pdf(d::PointMass, x::Real) = x == d.μ ? 1.0 : 0.0
Distributions.logpdf(d::PointMass, x::Real) = x == d.μ ? 0.0 : -Inf
Distributions.mean(d::PointMass) = d.μ
Distributions.std(d::PointMass) = 0
Distributions.var(d::PointMass) = 0

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

# Given α,β, get μ and σ
function betaMoments(dist::Distributions.Beta)
    α = dist.α
    β = dist.β

    μ = α / (α + β)
    σ = sqrt( α*β/((α+β)^2 * (α+β+1)) )
    return μ, σ
end

function gammaMoments(dist::Distributions.Gamma)
    α = dist.α
    θ = dist.θ

    μ = α / θ
    σ = α / θ^2
    return μ, σ
end


function Distributions.pdf(d::RootInverseGamma, x::Real)
    (ν, τ) = params(d)
    return 2 * (ν*τ^2/2)^(ν/2) * exp((-ν*τ^2)/(2x^2)) / gamma(ν/2) / x^(ν+1)
end

function Distributions.logpdf(d::RootInverseGamma, x::Real)
    (ν, τ) = params(d)
    return log(2) - log(gamma(ν/2)) + (ν/2)*log(ν*τ^2/2) - ((ν+1)/2)*log(x^2) - ν*τ^2/(2x^2)
end



# Degenerate multivariate normal
# en.wikipedia.org/wiki/Multivariate_normal_distribution#Degenerate_case
type DegenerateMvNormal <: Distribution{Multivariate, Continuous}
    μ::Vector          # mean
    σ::Matrix          # standard deviation
    rank::Int64        # rank
end

function DegenerateMvNormal(μ::Vector,σ::Matrix)
    return DegenerateMvNormal(μ,σ,rank(σ))
end

Base.length(d::DegenerateMvNormal) = length(d.μ)

# Generate a draw from d with variance optionally scaled by cc^2
function Distributions.rand{T<:AbstractFloat}(d::DegenerateMvNormal; cc::T = 1.0)
    return d.μ + cc*d.σ*randn(length(d))
end
