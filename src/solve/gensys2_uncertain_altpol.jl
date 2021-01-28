"""
```
gensys2_uncertain_altpol(prob_vec::AbstractVector{S}, Th::AbstractMatrix{S}, Ch::AbstractVector{S},
                     Tzlbs::Vector{Matrix{S}}, Rzlbs::Vector{Matrix{S}}, Czlbs::Vector{Vector{S}},
                     Γ0_til::AbstractMatrix{S}, Γ1_til::AbstractMatrix{S}, Γ2_til::AbstractMatrix{S},
                     C_til::AbstractVector{S}, Ψ_til::AbstractMatrix{S}) where {S <: Real}

gensys2_uncertain_altpol(prob_vec::AbstractVector{S}, Th::Vector{Matrix{Float64}}, Ch::Vector{Vector{Float64}},
                     Tzlbs::Vector{Matrix{S}}, Rzlbs::Vector{Matrix{S}}, Czlbs::Vector{Vector{S}},
                     Γ0_til::AbstractMatrix{S}, Γ1_til::AbstractMatrix{S}, Γ2_til::AbstractMatrix{S},
                     C_til::AbstractVector{S}, Ψ_til::AbstractMatrix{S}) where {S <: Real}
```

calculates the transition matrices when there is uncertainty about whether a temporary altpolicy. In particular,
- the `Th` and `Ch` are the "historical" matrices associated with policy the central bank
  no longer wants to use (but agents believe may still be occurring); in the case Th is a vector of matrices
  and Ch a vector of vectors, they are the matrices associated with all policies except the
  alternative_policy that agents believe may occur.
- the `til` matrices are the predictable form of the gensys matrices under the
  policy the central bank claims to implement;
- `T_alts`, `R_alts`, and `C_alts` are the gensys2 matrices with a temporary altpolicy imposed over
  a pre-specified horizon with the belief that the claimed policy actually occurs after the temporary altpolicy  ends.

As an example, the Federal Reserve in 2020-Q3 switched to flexible AIT. The "historical" policy is
a Taylor-style rule, and the `til` matrices correspond to the predictable form of flexible AIT.
The `T_alts`, `R_alts`, and `C_alts` include temporary switches to a ZLB specified by the `regime_eqcond_info`.
"""
function gensys2_uncertain_altpol(prob_vec::AbstractVector{S}, Th::Matrix{S}, Ch::Vector{S},
                              T_alts::Vector{Matrix{S}}, R_alts::Vector{Matrix{S}}, C_alts::Vector{Vector{S}},
                              Γ0_til::AbstractMatrix{S}, Γ1_til::AbstractMatrix{S}, Γ2_til::AbstractMatrix{S},
                              C_til::AbstractVector{S}, Ψ_til::AbstractMatrix{S}) where {S <: Real}

    Tbars = Vector{Matrix{S}}(undef, length(T_alts) + 1)
    Cbars = Vector{Vector{S}}(undef, length(T_alts) + 1)
    Tout  = Vector{Matrix{S}}(undef, length(T_alts) + 1)
    Rout  = Vector{Matrix{S}}(undef, length(T_alts) + 1)
    Cout  = Vector{Vector{S}}(undef, length(T_alts) + 1)

    # Calculate "uncertain" temporary altpolicy matrices and back out the transition equation
    for i in 1:length(T_alts)
        Tbars[i] = prob_vec[1] * T_alts[i] + prob_vec[2] * Th # it is assumed T_alts is a vector of the T_{t + 1}^{(temporary altpolicy)} matrices
        Cbars[i] = prob_vec[1] * C_alts[i] + prob_vec[2] * Ch

        Tout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ Γ1_til
        Rout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ Ψ_til
        Cout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ (C_til  - Γ2_til * Cbars[i])
    end

    # Add terminal condition
    Tout[end] = T_alts[end]
    Rout[end] = R_alts[end]
    Cout[end] = C_alts[end]

    return Tout, Rout, Cout
end

# With time-varying credibility
function gensys2_uncertain_altpol(prob_vec::Vector{Vector{S}}, Th::AbstractMatrix{S}, Ch::AbstractVector{S},
                              T_alts::Vector{Matrix{S}}, R_alts::Vector{Matrix{S}}, C_alts::Vector{Vector{S}},
                              Γ0_til::AbstractMatrix{S}, Γ1_til::AbstractMatrix{S}, Γ2_til::AbstractMatrix{S},
                              C_til::AbstractVector{S}, Ψ_til::AbstractMatrix{S}) where {S <: Real}

    Tbars = Vector{Matrix{S}}(undef, length(T_alts) + 1)
    Cbars = Vector{Vector{S}}(undef, length(T_alts) + 1)
    Tout  = Vector{Matrix{S}}(undef, length(T_alts) + 1)
    Rout  = Vector{Matrix{S}}(undef, length(T_alts) + 1)
    Cout  = Vector{Vector{S}}(undef, length(T_alts) + 1)

    # Calculate "uncertain" temporary altpolicy matrices and back out the transition equation
    for i in 1:length(T_alts)
        Tbars[i] = prob_vec[i][1] * T_alts[i] + (1.0-prob_vec[i][1]) * Th # it is assumed T_alts is a vector of the T_{t + 1}^{(temporary altpolicy)} matrices
        Cbars[i] = prob_vec[i][1] * C_alts[i] + (1.0-prob_vec[i][1]) * Ch

        Tout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ Γ1_til
        Rout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ Ψ_til
        Cout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ (C_til  - Γ2_til * Cbars[i])
    end

    # Add terminal condition
    Tout[end] = T_alts[end]
    Rout[end] = R_alts[end]
    Cout[end] = C_alts[end]

    return Tout, Rout, Cout
end

# Same 2 functions as above but where Th and Ch are vectors of matrices/vectors.

function gensys2_uncertain_altpol(prob_vec::AbstractVector{S},
                              Th::Vector{Matrix{Float64}}, Ch::Vector{Vector{Float64}},
                              T_alts::Vector{Matrix{S}}, R_alts::Vector{Matrix{S}}, C_alts::Vector{Vector{S}},
                              Γ0_til::AbstractMatrix{S}, Γ1_til::AbstractMatrix{S}, Γ2_til::AbstractMatrix{S},
                              C_til::AbstractVector{S}, Ψ_til::AbstractMatrix{S}) where {S <: Real}

    Tbars = Vector{Matrix{S}}(undef, length(T_alts) + 1)
    Cbars = Vector{Vector{S}}(undef, length(T_alts) + 1)
    Tout  = Vector{Matrix{S}}(undef, length(T_alts) + 1)
    Rout  = Vector{Matrix{S}}(undef, length(T_alts) + 1)
    Cout  = Vector{Vector{S}}(undef, length(T_alts) + 1)

    # Calculate "uncertain" temporary altpolicy matrices and back out the transition equation
    for i in 1:length(T_alts)
        Tbars[i] = prob_vec[1] * T_alts[i] .+ sum([prob_vec[j+1] .* Th[j] for j in 1:length(Th)]) # it is assumed T_alts is a vector of the T_{t + 1}^{(temporary altpolicy)} matrices
        Cbars[i] = prob_vec[1] * C_alts[i] .+ sum([prob_vec[j+1] .* Ch[j] for j in 1:length(Ch)])

        Tout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ Γ1_til
        Rout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ Ψ_til
        Cout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ (C_til  - Γ2_til * Cbars[i])
    end

    # Add terminal condition
    Tout[end] = T_alts[end]
    Rout[end] = R_alts[end]
    Cout[end] = C_alts[end]

    return Tout, Rout, Cout
end

# With time-varying credibility
function gensys2_uncertain_altpol(prob_vec::Vector{Vector{S}},Th::Vector{Matrix{Float64}},Ch::Vector{Vector{Float64}},
                              T_alts::Vector{Matrix{S}}, R_alts::Vector{Matrix{S}}, C_alts::Vector{Vector{S}},
                              Γ0_til::AbstractMatrix{S}, Γ1_til::AbstractMatrix{S}, Γ2_til::AbstractMatrix{S},
                              C_til::AbstractVector{S}, Ψ_til::AbstractMatrix{S}) where {S <: Real}

    Tbars = Vector{Matrix{S}}(undef, length(T_alts) + 1)
    Cbars = Vector{Vector{S}}(undef, length(T_alts) + 1)
    Tout  = Vector{Matrix{S}}(undef, length(T_alts) + 1)
    Rout  = Vector{Matrix{S}}(undef, length(T_alts) + 1)
    Cout  = Vector{Vector{S}}(undef, length(T_alts) + 1)

    # Calculate "uncertain" temporary altpolicy matrices and back out the transition equation
    for i in 1:length(T_alts)
        Tbars[i] = prob_vec[i][1] .* T_alts[i] .+ sum([prob_vec[i][j+1] .* Th[j] for j in 1:length(Th)]) # it is assumed T_alts is a vector of the T_{t + 1}^{(temporary altpolicy)} matrices
        Cbars[i] = prob_vec[i][1] * C_alts[i] .+ sum([prob_vec[i][j+1] .* Ch[j] for j in 1:length(Ch)])

        Tout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ Γ1_til
        Rout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ Ψ_til
        Cout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ (C_til  - Γ2_til * Cbars[i])
    end

    # Add terminal condition
    Tout[end] = T_alts[end]
    Rout[end] = R_alts[end]
    Cout[end] = C_alts[end]

    return Tout, Rout, Cout
end
