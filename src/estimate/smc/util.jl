"""
```
mvnormal_mixture_draw{T<:AbstractFloat}(p, Σ; cc, α, d_prop)
```

Create a `DegenerateMvNormal` distribution object, `d`, from a parameter vector, `p`, and a
covariance matrix, `Σ`.

Generate a draw from the mixture distribution of the `DegenerateMvNormal` scaled by `cc^2`
and with mixture proportion `α`, a `DegenerateMvNormal` centered at the same mean, but with a
covariance matrix of the diagonal entries of `Σ` scaled by `cc^2` with mixture
proportion `(1 - α)/2`, and an additional proposed distribution with the same covariance
matrix as `d` but centered at the new proposed mean, `p_prop`, scaled by `cc^2`, and with mixture proportion `(1 - α)/2`.

If no `p_prop` is given, but an `α` is specified, then the mixture will consist of `α` of
the standard distribution and `(1 - α)` of the diagonalized covariance distribution.

### Arguments
`p`: The mean of the desired distribution
`Σ`: The standard deviation of the desired distribution

"""
function mvnormal_mixture_draw{T<:AbstractFloat}(p::Vector{T}, Σ::Matrix{T};
                                                 cc::T = 1.0, α::T = 1.,
                                                 p_prop::Vector{T} = zeros(length(p)))
    @assert 0 <= α <= 1
    d = DegenerateMvNormal(p, Σ)
    d_diag = DegenerateMvNormal(p, diagm(diag(Σ)))
    d_prop = p_prop == zeros(length(p)) ? d_diag : DegenerateMvNormal(p_prop, Σ)

    normal_component = α*(d.μ + cc*d.σ*randn(length(d)))
    diag_component   = (1 - α)/2*(d_diag.μ + cc*d_diag.σ*randn(length(d_diag)))
    proposal_component   = (1 - α)/2*(d_prop.μ + cc*d_prop.σ*randn(length(d_prop)))

    return normal_component + diag_component + proposal_component
end
