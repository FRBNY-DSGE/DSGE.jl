using Distributions


# Distribution of fixed parameters
type PointMass <: Distribution{Univariate, Continuous}
    μ::Float64
end

function Distributions.pdf(d::PointMass, x::Float64)
    x == d.μ ? 1.0 : 0.0
end



# Given mean and SD, calculate α and β, return distribution
function Beta(μ::Float64, σ::Float64)
    α = (1-μ) * μ^2 / σ^2 - μ
    β = α * (1/μ - 1)
    return Distributions.Beta(α, β)
end

function Gamma(μ::Float64, σ::Float64)
    β = σ^2 / μ
    α = μ / β
    return Distributions.Gamma(α, β)
end

# ν and σ² parameterize the scaled inverse chi-squared distribution, a different parameterization of the inverse gamma distribution
function InverseGamma(σ::Float64, ν::Float64)
    α = ν / 2
    β = ν * σ^2 / 2
    return Distributions.InverseGamma(α, β)
end    
