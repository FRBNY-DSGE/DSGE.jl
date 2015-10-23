#=
"""
This file defines additional functions to return objects of type Distribution. This is
necessary because we specify prior distributions wrt mean and SD
(for beta and gamma-distributed parameters) and ν and σ (for inverse gamma-distributed
parameters). Note these functions are NOT new methods for the Distributions.Beta, etc.
functions, but rather new functions with the same names.
"""
=#

using Distributions, Compat, DSGE
import Distributions: params, mean, std, pdf, logpdf, rand
import Base: length

export PointMass, Beta, Gamma, RootInverseGamma, DegenerateMvNormal

#=
doc"""
PointMass <: Distribution{Univariate, Continuous}

Fixed parameters are drawn from a PointMass distribution, and thus always take on the same value.
"""
=#
type PointMass <: Distribution{Univariate, Continuous}
    μ::Float64
end

Distributions.pdf(d::PointMass, x::Real) = x == d.μ ? 1.0 : 0.0
Distributions.logpdf(d::PointMass, x::Real) = x == d.μ ? 0.0 : -Inf
Distributions.mean(d::PointMass) = d.μ
Distributions.std(d::PointMass) = 0
Distributions.var(d::PointMass) = 0

#=
doc"""
BetaAlt(μ::Real, σ::Real)

### Parameters
`μ::Real`: The mean of the desired distribution 
`σ::Real`: The standard deviation of the desired distribution 

### Description:
Given μ and σ, calculate α and β and return a Distributions.Beta Distribution object.
"""
=#
function BetaAlt(μ::Real, σ::Real)
    α = (1-μ) * μ^2 / σ^2 - μ
    β = α * (1/μ - 1)
    return Distributions.Beta(α, β)
end


#=
doc"""
GammaAlt(μ::Real, σ::Real)

### Parameters
`μ::Real`: The mean of the desired distribution 
`σ::Real`: The standard deviation of the desired distribution 

### Description:
Given μ and σ, calculate α and β and return a Distributions.Gamma object.
"""
=#
function GammaAlt(μ::Real, σ::Real)
    β = σ^2 / μ
    α = μ / β
    return Distributions.Gamma(α, β)
end

#=
doc"""
type RootInverseGamma <: Distribution{Univariate, Continuous}

If x ~ RootInverseGamma(ν, τ²), then
  x² ~ ScaledInverseChiSquared(ν, τ²)
  x² ~ InverseGamma(ν/2, ντ²/2)


"""
=#
type RootInverseGamma <: Distribution{Univariate, Continuous}
    ν::Float64
    τ::Float64
end

Distributions.params(d::RootInverseGamma) = (d.ν, d.τ)

#=
doc"""
betaMoments(dist::Distributions.Beta)

Compute the mean μ and standard deviation σ of a Distributions.Beta object. 
"""
=#
function betaMoments(dist::Distributions.Beta)
    α = dist.α
    β = dist.β

    μ = α / (α + β)
    σ = sqrt( α*β/((α+β)^2 * (α+β+1)) )
    return μ, σ
end

#=
doc"""
gammaMoments(dist::Distributions.Gamma)

Compute the mean μ and standard deviation σ of a Distributions.Gamma object. 
"""
=#
function gammaMoments(dist::Distributions.Gamma)
    α = dist.α
    θ = dist.θ

    μ = α / θ
    σ = α / θ^2
    return μ, σ
end

#=
doc"""
Distributions.pdf(d::RootInverseGamma, x::Real)

Compute the pdf of a RootInverseGamma distribution at x.
"""
=#
function Distributions.pdf(d::RootInverseGamma, x::Real)
    (ν, τ) = params(d)
    return 2 * (ν*τ^2/2)^(ν/2) * exp((-ν*τ^2)/(2x^2)) / gamma(ν/2) / x^(ν+1)
end

#=
doc"""
Distributions.logpdf(d::RootInverseGamma, x::Real)

Compute the log pdf of a RootInverseGamma distribution at x.
"""
=#
function Distributions.logpdf(d::RootInverseGamma, x::Real)
    (ν, τ) = params(d)
    return log(2) - log(gamma(ν/2)) + (ν/2)*log(ν*τ^2/2) - ((ν+1)/2)*log(x^2) - ν*τ^2/(2x^2)
end


#=
doc"""
DegenerateMvNormal <: Distribution{Multivariate, Continuous}

The DegenerateMvNormal type implements a degenerate multivariate normal distribution as a subtype of Distribution.
en.wikipedia.org/wiki/Multivariate_normal_distribution#Degenerate_case
"""
=#
type DegenerateMvNormal <: Distribution{Multivariate, Continuous}
    μ::Vector          # mean
    σ::Matrix          # standard deviation
    rank::Int64        # rank
end

#=
doc"""
DegenerateMvNormal(μ::Vector,σ::Matrix)

Returns a DegenerateMvNormal type with mean vector μ and covariance matrix σ
"""
=#
function DegenerateMvNormal(μ::Vector,σ::Matrix)
    return DegenerateMvNormal(μ,σ,rank(σ))
end

#=
doc"""
length(d::DegenerateMvNormal)

Returns the dimension of d.
"""
=#
Base.length(d::DegenerateMvNormal) = length(d.μ)

#=
doc"""
Distributions.rand{T<:AbstractFloat}(d::DegenerateMvNormal; cc::T = 1.0)

Generate a draw from d with variance optionally scaled by cc^2.
"""
=#
function Distributions.rand{T<:AbstractFloat}(d::DegenerateMvNormal; cc::T = 1.0)
    return d.μ + cc*d.σ*randn(length(d))
end

## #=
## doc"""
## Distributions.rand{T<:AbstractFloat}(d::DegenerateMvNormal; cc::T = 1.0)

## Generate a draw from d with variance optionally scaled by cc^2.
## """
## =#
## function Distributions.rand{T<:AbstractFloat, U<:AbstractDSGEModel}(d::DegenerateMvNormal, m::U; cc::T = 1.0)
##     return d.μ + cc*d.σ*randn(m.rng, length(d))
## end

