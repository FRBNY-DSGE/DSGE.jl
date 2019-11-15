"""
```
transition(m::PoolModel{T}) where {T<:AbstractFloat}
```

Assign transition equation

```
x_t = (1-ρ) μ + ρ x_{t-1} + sqrt{1 - ρ^2} σ ϵ_t
ϵ_t ∼ iid N(0,1), x_0 ∼ N(μ,σ^2)
λ_t = Φ(x_t)
```
where Φ(⋅) is the cdf of a N(0,1) random variable, F_ϵ is the distribution of ϵ_t,
and F_λ is the distribution of λ(x_0).
"""
function transition(m::PoolModel{T}) where {T<:AbstractFloat}
    @inline Φ(x::Vector{Float64}, ϵ::Vector{Float64}) = abs.([0;1] .-
                                                             (cdf.(Normal(), (1 - m[:ρ].value) .*
                                                                   m[:μ].value .+ m[:ρ].value .*
                                                                   quantile(Normal(),x[1]) .+
                                                                   sqrt(1 - m[:ρ].value^2) .*
                                                                   m[:σ].value .* ϵ)))
    F_ϵ = Normal(0.,1.)
    F_λ = Uniform(0.,1.)
    return Φ, F_ϵ, F_λ
end
