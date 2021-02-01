"""
```
gensys2_uncertain_altpol(prob_vec::AbstractVector{S}, T_alt::AbstractMatrix{S}, C_alt::AbstractVector{S},
                     T_impl::Vector{Matrix{S}}, R_impl::Vector{Matrix{S}}, C_impl::Vector{Vector{S}},
                     Γ0_til::AbstractMatrix{S}, Γ1_til::AbstractMatrix{S}, Γ2_til::AbstractMatrix{S},
                     C_til::AbstractVector{S}, Ψ_til::AbstractMatrix{S}) where {S <: Real}

gensys2_uncertain_altpol(prob_vec::AbstractVector{S}, T_alt::Vector{Matrix{S}}, C_alt::Vector{Vector{S}},
                     T_impl::Vector{Matrix{S}}, R_impl::Vector{Matrix{S}}, C_impl::Vector{Vector{S}},
                     Γ0_til::AbstractMatrix{S}, Γ1_til::AbstractMatrix{S}, Γ2_til::AbstractMatrix{S},
                     C_til::AbstractVector{S}, Ψ_til::AbstractMatrix{S}) where {S <: Real}
```

calculates the transition matrices when there is a temporary alternative policy with imperfect awareness.
In particular,
- the `T_alt` and `C_alt` are the matrices associated with the other alternative policies
    agents believe may occur. If T_alt is a vector of matrices and C_alt a vector of vectors,
    then there are multiple alternative policies agents believe may occur.
- the `til` matrices are the predictable form of the gensys matrices under the
    policy the central bank actually implements during the horizon of the temporary alternative policy.
- `T_impl`, `R_impl`, and `C_impl` are the gensys2 matrices of the temporary alternative policy
    that is actually implemented.

As an example, the Federal Reserve in 2020-Q3 switched to flexible AIT. T_alte "historical" policy is
a Taylor-style rule, and the `til` matrices correspond to the predictable form of flexible AIT.
T_alte `T_impl`, `R_impl`, and `C_impl` include temporary switches to a ZLB specified by the `regime_eqcond_info`.

Note that it is still assumed that only one policy occurs during the temporary alternative policy horizon.
For example, we cannot have a temporary ZLB, followed by temporary AIT, before switching back to a Taylor Rule.
To generalize further, we need to make `Γ0_til`, etc., a vector of matrices/vectors and figure out a way
to easily pass these conditions in.
"""
function gensys2_uncertain_altpol(prob_vec::AbstractVector{S}, T_alt::Matrix{S}, C_alt::Vector{S},
                                  T_impl::Vector{Matrix{S}}, R_impl::Vector{Matrix{S}}, C_impl::Vector{Vector{S}},
                                  Γ0_til::AbstractMatrix{S}, Γ1_til::AbstractMatrix{S}, Γ2_til::AbstractMatrix{S},
                                  C_til::AbstractVector{S}, Ψ_til::AbstractMatrix{S}) where {S <: Real}

    Tbars = Vector{Matrix{S}}(undef, length(T_impl) + 1)
    Cbars = Vector{Vector{S}}(undef, length(T_impl) + 1)
    Tout  = Vector{Matrix{S}}(undef, length(T_impl) + 1)
    Rout  = Vector{Matrix{S}}(undef, length(T_impl) + 1)
    Cout  = Vector{Vector{S}}(undef, length(T_impl) + 1)

    # Calculate "uncertain" temporary altpolicy matrices and back out the transition equation
    for i in 1:length(T_impl)
        Tbars[i] = prob_vec[1] * T_impl[i] + prob_vec[2] * T_alt # it is assumed T_impl is a vector of the
        Cbars[i] = prob_vec[1] * C_impl[i] + prob_vec[2] * C_alt # T_{t + 1}^{(temporary altpolicy)} matrices
        # hence Tbars will be a vector of T_t matrices
        Tout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ Γ1_til
        Rout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ Ψ_til
        Cout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ (C_til  - Γ2_til * Cbars[i])
    end

    # Add terminal condition
    Tout[end] = T_impl[end]
    Rout[end] = R_impl[end]
    Cout[end] = C_impl[end]

    return Tout, Rout, Cout
end

# With time-varying credibility
function gensys2_uncertain_altpol(prob_vec::Vector{Vector{S}}, T_alt::AbstractMatrix{S}, C_alt::AbstractVector{S},
                                  T_impl::Vector{Matrix{S}}, R_impl::Vector{Matrix{S}}, C_impl::Vector{Vector{S}},
                                  Γ0_til::AbstractMatrix{S}, Γ1_til::AbstractMatrix{S}, Γ2_til::AbstractMatrix{S},
                                  C_til::AbstractVector{S}, Ψ_til::AbstractMatrix{S}) where {S <: Real}

    Tbars = Vector{Matrix{S}}(undef, length(T_impl) + 1)
    Cbars = Vector{Vector{S}}(undef, length(T_impl) + 1)
    Tout  = Vector{Matrix{S}}(undef, length(T_impl) + 1)
    Rout  = Vector{Matrix{S}}(undef, length(T_impl) + 1)
    Cout  = Vector{Vector{S}}(undef, length(T_impl) + 1)

    # Calculate "uncertain" temporary altpolicy matrices and back out the transition equation
    for i in 1:length(T_impl)
        Tbars[i] = prob_vec[i][1] * T_impl[i] + (1.0-prob_vec[i][1]) * T_alt # it is assumed T_impl is a vector of the
        Cbars[i] = prob_vec[i][1] * C_impl[i] + (1.0-prob_vec[i][1]) * C_alt # T_{t + 1}^{(temporary altpolicy)} matrices
        # hence Tbars will be a vector of T_t matrices
        Tout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ Γ1_til
        Rout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ Ψ_til
        Cout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ (C_til  - Γ2_til * Cbars[i])
    end

    # Add terminal condition
    Tout[end] = T_impl[end]
    Rout[end] = R_impl[end]
    Cout[end] = C_impl[end]

    return Tout, Rout, Cout
end

# Same 2 functions as above but where T_alt and C_alt are vectors of matrices/vectors.
function gensys2_uncertain_altpol(prob_vec::AbstractVector{S},
                                  T_alt::Vector{Matrix{S}}, C_alt::Vector{Vector{S}},
                                  T_impl::Vector{Matrix{S}}, R_impl::Vector{Matrix{S}}, C_impl::Vector{Vector{S}},
                                  Γ0_til::AbstractMatrix{S}, Γ1_til::AbstractMatrix{S}, Γ2_til::AbstractMatrix{S},
                                  C_til::AbstractVector{S}, Ψ_til::AbstractMatrix{S}) where {S <: Real}

    Tbars = Vector{Matrix{S}}(undef, length(T_impl) + 1)
    Cbars = Vector{Vector{S}}(undef, length(T_impl) + 1)
    Tout  = Vector{Matrix{S}}(undef, length(T_impl) + 1)
    Rout  = Vector{Matrix{S}}(undef, length(T_impl) + 1)
    Cout  = Vector{Vector{S}}(undef, length(T_impl) + 1)

    # Calculate "uncertain" temporary altpolicy matrices and back out the transition equation
    for i in 1:length(T_impl)
        Tbars[i] = prob_vec[1] * T_impl[i] .+ sum([prob_vec[j+1] .* T_alt[j] for j in 1:length(T_alt)]) # it is assumed T_impl is a vector of the
        Cbars[i] = prob_vec[1] * C_impl[i] .+ sum([prob_vec[j+1] .* C_alt[j] for j in 1:length(C_alt)]) # T_{t + 1}^{(temporary altpolicy)} matrices

        Tout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ Γ1_til
        Rout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ Ψ_til
        Cout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ (C_til  - Γ2_til * Cbars[i])
    end

    # Add terminal condition
    Tout[end] = T_impl[end]
    Rout[end] = R_impl[end]
    Cout[end] = C_impl[end]

    return Tout, Rout, Cout
end

# With time-varying credibility
function gensys2_uncertain_altpol(prob_vec::Vector{Vector{S}},T_alt::Vector{Matrix{S}},C_alt::Vector{Vector{S}},
                                  T_impl::Vector{Matrix{S}}, R_impl::Vector{Matrix{S}}, C_impl::Vector{Vector{S}},
                                  Γ0_til::AbstractMatrix{S}, Γ1_til::AbstractMatrix{S}, Γ2_til::AbstractMatrix{S},
                                  C_til::AbstractVector{S}, Ψ_til::AbstractMatrix{S}) where {S <: Real}

    Tbars = Vector{Matrix{S}}(undef, length(T_impl) + 1)
    Cbars = Vector{Vector{S}}(undef, length(T_impl) + 1)
    Tout  = Vector{Matrix{S}}(undef, length(T_impl) + 1)
    Rout  = Vector{Matrix{S}}(undef, length(T_impl) + 1)
    Cout  = Vector{Vector{S}}(undef, length(T_impl) + 1)

    # Calculate "uncertain" temporary altpolicy matrices and back out the transition equation
    for i in 1:length(T_impl)
        Tbars[i] = prob_vec[i][1] .* T_impl[i] .+ sum([prob_vec[i][j+1] .* T_alt[j] for j in 1:length(T_alt)]) # it is assumed T_impl is a vector of the
        Cbars[i] = prob_vec[i][1] * C_impl[i] .+ sum([prob_vec[i][j+1] .* C_alt[j] for j in 1:length(C_alt)])  # T_{t + 1}^{(temporary altpolicy)} matrices

        Tout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ Γ1_til
        Rout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ Ψ_til
        Cout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ (C_til  - Γ2_til * Cbars[i])
    end

    # Add terminal condition
    Tout[end] = T_impl[end]
    Rout[end] = R_impl[end]
    Cout[end] = C_impl[end]

    return Tout, Rout, Cout
end

# Same 2 functions as above but where T_alt and C_alt are
# are a Vector{Union{Matrix/Vector{S}, Vector{Matrix/Vector{S}}}} (to allow heterogeneous Vector and Matrix types)
function gensys2_uncertain_altpol(prob_vec::AbstractVector{S}, gensys2_regimes::UnitRange{Int64},
                                  T_alt::Vector{Union{Matrix{S}, Vector{Matrix{S}}}},
                                  C_alt::Vector{Union{Vector{S}, Vector{Vector{S}}}}, is_altpol::Vector{Bool},
                                  T_impl::Vector{Matrix{S}}, R_impl::Vector{Matrix{S}}, C_impl::Vector{Vector{S}},
                                  Γ0_til::AbstractMatrix{S}, Γ1_til::AbstractMatrix{S}, Γ2_til::AbstractMatrix{S},
                                  C_til::AbstractVector{S}, Ψ_til::AbstractMatrix{S}) where {S <: Real}

    Tbars = Vector{Matrix{S}}(undef, length(T_impl) + 1)
    Cbars = Vector{Vector{S}}(undef, length(T_impl) + 1)
    Tout  = Vector{Matrix{S}}(undef, length(T_impl) + 1)
    Rout  = Vector{Matrix{S}}(undef, length(T_impl) + 1)
    Cout  = Vector{Vector{S}}(undef, length(T_impl) + 1)

    # Calculate "uncertain" temporary altpolicy matrices and back out the transition equation
    n_alt = length(prob_vec) - 1
    for (i, reg) in zip(1:length(T_impl), gensys2_regimes[1:end - 1])
        # it is assumed T_impl is a vector of the T_{t + 1}^{(temporary altpolicy)} matrices,
        # hence Tbars will be a vector of T_t matrices
        Tbars[i] = prob_vec[1] * T_impl[i] .+
        sum([prob_vec[j+1] .* (is_altpol[j] ? T_alt[j] : T_alt[j][reg])
             for j in 1:n_alt])
        Cbars[i] = prob_vec[1] * C_impl[i] .+
        sum([prob_vec[j+1] .* (is_altpol[j] ? C_alt[j] : C_alt[j][reg])
             for j in 1:n_alt])

        Tout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ Γ1_til
        Rout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ Ψ_til
        Cout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ (C_til  - Γ2_til * Cbars[i])
    end

    # Add terminal condition
    Tout[end] = T_impl[end]
    Rout[end] = R_impl[end]
    Cout[end] = C_impl[end]

    return Tout, Rout, Cout
end

# With time-varying credibility
function gensys2_uncertain_altpol(prob_vec::Vector{Vector{S}}, gensys2_regimes::UnitRange{Int64},
                                  T_alt::Vector{Union{Matrix{S}, Vector{Matrix{S}}}},
                                  C_alt::Vector{Union{Vector{S}, Vector{Vector{S}}}}, is_altpol::Vector{Bool},
                                  Γ0_til::AbstractMatrix{S}, Γ1_til::AbstractMatrix{S}, Γ2_til::AbstractMatrix{S},
                                  C_til::AbstractVector{S}, Ψ_til::AbstractMatrix{S}) where {S <: Real}

    Tbars = Vector{Matrix{S}}(undef, length(gensys2_regimes))
    Cbars = Vector{Vector{S}}(undef, length(gensys2_regimes))
    Tout  = Vector{Matrix{S}}(undef, length(gensys2_regimes))
    Rout  = Vector{Matrix{S}}(undef, length(gensys2_regimes))
    Cout  = Vector{Vector{S}}(undef, length(gensys2_regimes))

    # Calculate "uncertain" temporary altpolicy matrices and back out the transition equation
    n_alt = length(prob_vec[1])
    for (i, reg) in enumerate(gensys2_regimes[1:end - 1])
        # it is assumed T_impl is a vector of the T_{t + 1}^{(temporary altpolicy)} matrices,
        # hence Tbars will be a vector of T_t matrices
        Tbars[i] = prob_vec[i][1] * T_alt[1][reg] .+
            sum([prob_vec[i][j+1] .* (is_altpol[j] ? T_alt[j] : T_alt[j][reg])
                 for j in 2:n_alt])
        Cbars[i] = prob_vec[i][1] * C_alt[1][reg] .+
            sum([prob_vec[i][j+1] .* (is_altpol[j] ? C_alt[j] : C_alt[j][reg])
                 for j in 2:n_alt])

        Tout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ Γ1_til
        Rout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ Ψ_til
        Cout[i]  = (Γ2_til * Tbars[i] + Γ0_til) \ (C_til  - Γ2_til * Cbars[i])
    end

    return Tout, Rout, Cout
end
