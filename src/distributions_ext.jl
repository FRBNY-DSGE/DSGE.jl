#=
This file defines additional functions to return objects of type Distribution. This is
necessary because we specify prior distributions wrt mean and SD
(for beta and gamma-distributed parameters) and ν and σ (for inverse gamma-distributed
parameters). Note these functions are NOT new methods for the Distributions.Beta, etc.
functions, but rather new functions with the same names.
=#

import Distributions: params, mean, std, pdf, logpdf, rand
import Base: length, rank

"""
```
BetaAlt(μ::AbstractFloat, σ::AbstractFloat)
```

Given μ and σ, calculate α and β and return a Distributions.Beta Distribution object.

### Arguments
`μ`: The mean of the desired distribution
`σ`: The standard deviation of the desired distribution
"""
function BetaAlt(μ::AbstractFloat, σ::AbstractFloat)
    α = (1-μ) * μ^2 / σ^2 - μ
    β = α * (1/μ - 1)
    return Distributions.Beta(α, β)
end


"""
```
GammaAlt(μ::AbstractFloat, σ::AbstractFloat)
```

Given μ and σ, calculate α and β and return a Distributions.Gamma object.

### Arguments
`μ`: The mean of the desired distribution
`σ`: The standard deviation of the desired distribution
"""
function GammaAlt(μ::AbstractFloat, σ::AbstractFloat)
    β = σ^2 / μ
    α = μ / β
    return Distributions.Gamma(α, β)
end

"""
```
type RootInverseGamma <: Distribution{Univariate, Continuous}
```

If x ~ RootInverseGamma(ν, τ²), then
  x² ~ ScaledInverseChiSquared(ν, τ²)
  x² ~ InverseGamma(ν/2, ντ²/2)

ν represents the degrees of freedom.
"""
type RootInverseGamma <: Distribution{Univariate, Continuous}
    ν::Float64
    τ::Float64
end

Distributions.params(d::RootInverseGamma) = (d.ν, d.τ)

"""
```
Distributions.pdf(d::RootInverseGamma, x::AbstractFloat)
```

Compute the pdf of a RootInverseGamma distribution at x.
"""
function Distributions.pdf(d::RootInverseGamma, x::AbstractFloat)
    (ν, τ) = params(d)
    return 2 * (ν*τ^2/2)^(ν/2) * exp((-ν*τ^2)/(2x^2)) / gamma(ν/2) / x^(ν+1)
end

"""
```
Distributions.logpdf(d::RootInverseGamma, x::AbstractFloat)
```

Compute the log pdf of a RootInverseGamma distribution at x.
"""
function Distributions.logpdf(d::RootInverseGamma, x::AbstractFloat)
    (ν, τ) = params(d)
    return log(2) - log(gamma(ν/2)) + (ν/2)*log(ν*τ^2/2) - ((ν+1)/2)*log(x^2) - ν*τ^2/(2x^2)
end

"""
```
Distributions.rand{T<:AbstractFloat}(d::RootInverseGamma; cc::T = 1.0)
```

Generate a draw from d with variance optionally scaled by cc^2 (for a RootInverseGamma)
"""
function Distributions.rand{T<:AbstractFloat}(d::RootInverseGamma; cc::T = 1.0)
    return sqrt(d.ν * d.τ^2 / sum(randn(round(Int,d.ν)).^2))
end

"""
```
DegenerateMvNormal <: Distribution{Multivariate, Continuous}
```

The `DegenerateMvNormal` type implements a degenerate multivariate normal
distribution. The covariance matrix may not be full rank (hence degenerate).

See [Multivariate normal distribution - Degenerate case](en.wikipedia.org/wiki/Multivariate_normal_distribution#Degenerate_case).
"""
type DegenerateMvNormal <: Distribution{Multivariate, Continuous}
    μ::Vector          # mean
    σ::Matrix          # standard deviation
end

"""
```
rank(d::DegenerateMvNormal)
```

Returns the rank of `d.σ`.
"""
function rank(d::DegenerateMvNormal)
    return rank(d.σ)
end

"""
```
length(d::DegenerateMvNormal)
```

Returns the dimension of `d`.
"""
Base.length(d::DegenerateMvNormal) = length(d.μ)

"""
```
Distributions.rand{T<:AbstractFloat}(d::DegenerateMvNormal; cc::T = 1.0)
```

Generate a draw from `d` with variance optionally scaled by `cc^2`.
"""
function Distributions.rand{T<:AbstractFloat}(d::DegenerateMvNormal; cc::T = 1.0)
    return d.μ + cc*d.σ*randn(length(d))
end

"""
```
Distributions.logpdf{T<:AbstractFloat}(d::DegenerateMvNormal, v::Vector)
```

Evaluate the logpdf of `d` at `v` subsetting out the positive-definite sub-matrix of the
covariance matrix of `d` and the corresponding indices in `v`.
"""
function Distributions.logpdf{T<:AbstractFloat}(d::DegenerateMvNormal, v::Vector{T})
    inds = find(!iszero, [d.σ[:, i] for i in 1:length(v)])
    cov_mat = d.σ[inds, inds]*d.σ[inds, inds]'
    d_alt = MvNormal(d.μ[inds], nearest_spd(cov_mat))
    v_alt = v[inds]
    return logpdf(d_alt, v_alt)
end

"""
```
moments(dist::Uniform)
```

Compute the mean μ and standard deviation σ of a DSGE.Uniform object.
"""
function moments(dist::Uniform)
    return ((dist.a + dist.b)/2) , ((dist.b - dist.a)^2/12)
end

"""
```
moments(dist::DegenerateMvNormal)
```

Compute the mean μ and standard deviation σ of a DSGE.DegenerateMvNormal object.
"""
function moments(dist::DegenerateMvNormal)
    return dist.μ, dist.σ
end

"""
```
moments(dist::RootInverseGamma)
```

Return the RootInverseGamma parameters τ and ν.
"""
function moments(dist::RootInverseGamma)
    ν, τ = params(dist)

    return τ, ν


    # Compute the mean μ and standard deviation σ of a DSGE.RootInverseGamma object.

    # A Root Inverse Gamma / Nagasaki Scaled Chi2 Distribution's mean and standard deviation
    # can be computed as follows:

    #     μ = √(β) * Γ(α - 0.5) / Γ(α)
    #     σ = √(β / (α - 1) - μ²)

    # where α = ν/2 and β = τ²*ν/2.

    # α = ν/2
    # β = τ^2*ν/2

    # μ = β^(.5) * gamma(α - 0.5) / gamma(α)
    # σ = (β / (α - 1) - μ^2)^(0.5)

    # return μ, σ
end


"""
```
moments(dist::Distributions.Beta)
```

Compute the mean μ and standard deviation σ of a Distributions.Beta object.
"""
function moments(dist::Distributions.Beta)
    α = dist.α
    β = dist.β

    μ = α / (α + β)
    σ = sqrt( α*β/((α+β)^2 * (α+β+1)) )
    return μ, σ
end

"""
```
moments(dist::Distributions.Gamma)
```

Compute the mean μ and standard deviation σ of a Distributions.Gamma object.
"""
function moments(dist::Distributions.Gamma)
    α = dist.α
    θ = dist.θ

    μ = α / θ
    σ = α / θ^2
    return μ, σ
end

"""
```
moments(dist::Distributions.Normal)
```

Return the mean and std dev or a Distributions.Normal object.
"""
function moments(dist::Distributions.Normal)
    params(dist)
end

"""
```
Distributions.rand(d::DegenerateMvNormal, n::Int)
```

Generate `n` draws from `d`. This returns a matrix of size `(length(d), n)`,
where each column is a sample.
"""
function Distributions.rand(d::DegenerateMvNormal, n::Int)
    return d.μ .+ d.σ*randn(length(d), n)
end

"""
```
DegenerateDiagMvTDist <: Distribution{Multivariate, Continuous}
```

The `DegenerateDiagMvTDist` type implements a degenerate multivariate Student's t
distribution, where the covariance matrix is diagonal. The covariance matrix may
not be full rank (hence degenerate).
"""
type DegenerateDiagMvTDist <: Distribution{Multivariate, Continuous}
    μ::Vector          # mean
    σ::Matrix          # standard deviation
    ν::Int             # degrees of freedom

    function DegenerateDiagMvTDist(μ::Vector, σ::Matrix, ν::Int)
        ν > 0       ? nothing : error("Degrees of freedom (ν) must be positive")
        isdiag(σ^2) ? nothing : error("Covariance matrix (σ^2) must be diagonal")
        return new(μ, σ, ν)
    end
end

"""
```
rank(d::DegenerateDiagMvTDist)
```

Returns the rank of `d.σ`.
"""
function rank(d::DegenerateDiagMvTDist)
    return rank(d.σ)
end

"""
```
length(d::DegenerateDiagMvTDist)
```

Returns the dimension of `d`.
"""
Base.length(d::DegenerateDiagMvTDist) = length(d.μ)

"""
```
Distributions.rand(d::DegenerateDiagMvTDist)
```

Generate a draw from `d`.
"""
function Distributions.rand(d::DegenerateDiagMvTDist)
    return d.μ + d.σ*rand(TDist(d.ν), length(d))
end

"""
```
Distributions.rand(d::DegenerateDiagMvTDist, n::Int)
```

Generate `n` draws from `d`. This returns a matrix of size `(length(d), n)`,
where each column is a sample.
"""
function Distributions.rand(d::DegenerateDiagMvTDist, n::Int)
    return d.μ .+ d.σ*rand(TDist(d.ν), length(d), n)
end

# Compute the mean μ and standard deviation σ of a DSGE.RootInverseGamma object.

# A Root Inverse Gamma / Nagasaki Scaled Chi2 Distribution's mean and standard deviation
# can be computed as follows:

#     μ = √(β) * Γ(α - 0.5) / Γ(α)
#     σ = √(β / (α - 1) - μ²)

# where α = ν/2 and β = τ²*ν/2.
function mean(dist::RootInverseGamma)
    α = dist.ν/2
    β = dist.τ^2 * dist.ν/2

    μ = β^(.5) * gamma(α - 0.5) / gamma(α)
    return μ
end

function std(dist::RootInverseGamma)
    α = dist.ν/2
    β = dist.τ^2 * dist.ν/2

    σ = (β / (α - 1) - μ^2)^(0.5)
    return σ
end