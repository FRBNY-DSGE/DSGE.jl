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

gensys2_uncertain_altpol(prob_vec::Vector{S}, T_alt::Vector{Matrix{S}}, C_alt::Vector{Vector{S}},
                         T_impl::Vector{Matrix{S}}, R_impl::Vector{Matrix{S}}, C_impl::Vector{Vector{S}},
                         Γ0_til::AbstractMatrix{S}, Γ1_til::AbstractMatrix{S}, Γ2_til::AbstractMatrix{S},
                         C_til::AbstractVector{S}, Ψ_til::AbstractMatrix{S}) where {S <: Real}

gensys2_uncertain_altpol(prob_vec::Vector{Vector{S}}, T_alt::Vector{Matrix{S}}, C_alt::Vector{Vector{S}},
                         T_impl::Vector{Matrix{S}}, R_impl::Vector{Matrix{S}}, C_impl::Vector{Vector{S}},
                         Γ0_til::AbstractMatrix{S}, Γ1_til::AbstractMatrix{S}, Γ2_til::AbstractMatrix{S},
                         C_til::AbstractVector{S}, Ψ_til::AbstractMatrix{S}) where {S <: Real}

gensys2_uncertain_altpol(prob_vec::Vector{S}, gensys2_regimes::Vector{Int64}, inds::UnitRange{Int64},
                         systems::Vector{Union{System, RegimeSwitchingSystem}}, is_altpol::Vector{Bool},
                         Γ0_til::AbstractMatrix{S}, Γ1_til::AbstractMatrix{S}, Γ2_til::AbstractMatrix{S},
                         C_til::AbstractVector{S}, Ψ_til::AbstractMatrix{S}) where {S <: Real}

gensys2_uncertain_altpol(prob_vec::Vector{Vector{S}}, gensys2_regimes::Vector{Int64}, inds::UnitRange{Int64},
                         systems::Vector{Union{System, RegimeSwitchingSystem}}, is_altpol::Vector{Bool},
                         Γ0_til::AbstractMatrix{S}, Γ1_til::AbstractMatrix{S}, Γ2_til::AbstractMatrix{S},
                         C_til::AbstractVector{S}, Ψ_til::AbstractMatrix{S}) where {S <: Real}
```

calculates the transition matrices when there is a temporary alternative policy with imperfect awareness.
In particular, for the first four methods
- the `T_alt` and `C_alt` are the matrices associated with the other alternative policies
    agents believe may occur. If T_alt is a vector of matrices and C_alt a vector of vectors,
    then there are multiple alternative policies agents believe may occur.
- the `til` matrices are the predictable form of the gensys matrices under the
    policy the central bank actually implements during the horizon of the temporary policy.
- `T_impl`, `R_impl`, and `C_impl` are the gensys2 matrices of the temporary policy
    that is actually implemented.

As an example, the Federal Reserve in 2020-Q3 switched to flexible AIT. The "alternative" policy is
a Taylor-style rule, and the `til` matrices correspond to the predictable form of flexible AIT.
The `T_impl`, `R_impl`, and `C_impl` include temporary switches to a ZLB specified by `regime_eqcond_info`
followed by a switch to a permanent flexible AIT policy.

Note that it is still assumed that only one policy type occurs during the temporary alternative policy horizon.
For example, we cannot have a temporary ZLB, followed by temporary AIT, before switching back to a Taylor Rule.
To generalize further, we need to make `Γ0_til`, etc., a vector of matrices/vectors and figure out a way
to easily pass these conditions in. However, we can have regime dependence in the temporary rule,
e.g. a ZLB rule that has different ZLB values depending on the regime.

The last two methods allow the user to specify temporary policies as alternative policies
which agents believe may occur but are not implemented.
- `gensys2_regimes`: specifies the correct regimes (and take the form [first_gensys2_regime - 1, first_gensys2_regime, ..., last_gensys2_regime],
    where first_gensys2_regime is the first regime in which a temporary policy applies and
    last_gensys2_regime is the lift-off regime, i.e. the regime in which a permanent policy now applies
    instead of a temporary one).
- `inds`: specifies the required endogenous states for solving the model (i.e. states which are not added by augment_states)
- `systems`: specifies the alternative policies (which may or may not include state space systems
    with temporary policies) agents believe may occur but are not implemented.
- `is_altpol`: a vector specifying which of the alternative policiesare AltPolicy types
    rather than MultiPeriodAltPolicy. Note this vector excludes the implemented policy, so it has the same length of `systems`.
"""
function gensys2_uncertain_altpol(prob_vec::AbstractVector{S}, T_alt::Matrix{S}, C_alt::Vector{S},
                                  T_impl::Vector{Matrix{S}}, R_impl::Vector{Matrix{S}}, C_impl::Vector{Vector{S}},
                                  Γ0_til::AbstractMatrix{S}, Γ1_til::AbstractMatrix{S}, Γ2_til::AbstractMatrix{S},
                                  C_til::AbstractVector{S}, Ψ_til::AbstractMatrix{S}) where {S <: Real}

    Tout  = Vector{Matrix{S}}(undef, length(T_impl) + 1)
    Rout  = Vector{Matrix{S}}(undef, length(T_impl) + 1)
    Cout  = Vector{Vector{S}}(undef, length(T_impl) + 1)

    # Calculate "uncertain" temporary altpolicy matrices and back out the transition equation
    for i in 1:length(T_impl)
        Tbars = prob_vec[1] * T_impl[i] + prob_vec[2] * T_alt # it is assumed T_impl is a vector of the
        Cbars = prob_vec[1] * C_impl[i] + prob_vec[2] * C_alt # T_{t + 1}^{(temporary altpolicy)} matrices

        # hence Tbars will be a vector of T_t matrices
        tmp     = Γ2_til * Tbars + Γ0_til
        Tout[i] = tmp \ Γ1_til
        Rout[i] = tmp \ Ψ_til
        Cout[i] = tmp \ (C_til  - Γ2_til * Cbars)
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

    Tout  = Vector{Matrix{S}}(undef, length(T_impl) + 1)
    Rout  = Vector{Matrix{S}}(undef, length(T_impl) + 1)
    Cout  = Vector{Vector{S}}(undef, length(T_impl) + 1)

    # Calculate "uncertain" temporary altpolicy matrices and back out the transition equation
    for i in 1:length(T_impl)
        Tbars = prob_vec[i][1] * T_impl[i] + (1.0-prob_vec[i][1]) * T_alt # it is assumed T_impl is a vector of the
        Cbars = prob_vec[i][1] * C_impl[i] + (1.0-prob_vec[i][1]) * C_alt # T_{t + 1}^{(temporary altpolicy)} matrices

        # hence Tbars will be a vector of T_t matrices
        tmp     = Γ2_til * Tbars + Γ0_til
        Tout[i] = tmp \ Γ1_til
        Rout[i] = tmp \ Ψ_til
        Cout[i] = tmp \ (C_til  - Γ2_til * Cbars)
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

    Tout  = Vector{Matrix{S}}(undef, length(T_impl) + 1)
    Rout  = Vector{Matrix{S}}(undef, length(T_impl) + 1)
    Cout  = Vector{Vector{S}}(undef, length(T_impl) + 1)

    # Calculate "uncertain" temporary altpolicy matrices and back out the transition equation
    for i in 1:length(T_impl)
        Tbars = prob_vec[1] * T_impl[i] .+ sum([prob_vec[j+1] .* T_alt[j] for j in 1:length(T_alt)]) # it is assumed T_impl is a vector of the
        Cbars = prob_vec[1] * C_impl[i] .+ sum([prob_vec[j+1] .* C_alt[j] for j in 1:length(C_alt)]) # T_{t + 1}^{(temporary altpolicy)} matrices

        tmp     = Γ2_til * Tbars + Γ0_til
        Tout[i] = tmp \ Γ1_til
        Rout[i] = tmp \ Ψ_til
        Cout[i] = tmp \ (C_til  - Γ2_til * Cbars)
    end

    # Add terminal condition
    Tout[end] = T_impl[end]
    Rout[end] = R_impl[end]
    Cout[end] = C_impl[end]

    return Tout, Rout, Cout
end

# With time-varying credibility
function gensys2_uncertain_altpol(prob_vec::Vector{Vector{S}}, T_alt::Vector{Matrix{S}}, C_alt::Vector{Vector{S}},
                                  T_impl::Vector{Matrix{S}}, R_impl::Vector{Matrix{S}}, C_impl::Vector{Vector{S}},
                                  Γ0_til::AbstractMatrix{S}, Γ1_til::AbstractMatrix{S}, Γ2_til::AbstractMatrix{S},
                                  C_til::AbstractVector{S}, Ψ_til::AbstractMatrix{S}) where {S <: Real}

    Tout  = Vector{Matrix{S}}(undef, length(T_impl) + 1)
    Rout  = Vector{Matrix{S}}(undef, length(T_impl) + 1)
    Cout  = Vector{Vector{S}}(undef, length(T_impl) + 1)

    # Calculate "uncertain" temporary altpolicy matrices and back out the transition equation
    for i in 1:length(T_impl)
        Tbars = prob_vec[i][1] .* T_impl[i] .+ sum([prob_vec[i][j+1] .* T_alt[j] for j in 1:length(T_alt)]) # it is assumed T_impl is a vector of the
        Cbars = prob_vec[i][1] * C_impl[i] .+ sum([prob_vec[i][j+1] .* C_alt[j] for j in 1:length(C_alt)])  # T_{t + 1}^{(temporary altpolicy)} matr

        tmp     = Γ2_til * Tbars + Γ0_til
        Tout[i] = tmp \ Γ1_til
        Rout[i] = tmp \ Ψ_til
        Cout[i] = tmp \ (C_til  - Γ2_til * Cbars)
    end

    # Add terminal condition
    Tout[end] = T_impl[end]
    Rout[end] = R_impl[end]
    Cout[end] = C_impl[end]

    return Tout, Rout, Cout
end

# Same 2 functions as above but where T_alt and C_alt replaced by
# systems, a vector of perfectly credible alternative policies, including the implemented one
function gensys2_uncertain_altpol(prob_vec::AbstractVector{S}, gensys2_regimes::UnitRange{Int64}, inds::UnitRange{Int64},
                                  systems::Vector{Union{System, RegimeSwitchingSystem}}, is_altpol::Vector{Bool},
                                  T_impl::Vector{Matrix{S}}, R_impl::Vector{Matrix{S}}, C_impl::Vector{Vector{S}},
                                  Γ0_til::AbstractMatrix{S}, Γ1_til::AbstractMatrix{S}, Γ2_til::AbstractMatrix{S},
                                  C_til::AbstractVector{S}, Ψ_til::AbstractMatrix{S}) where {S <: Real}

    Tout  = Vector{Matrix{S}}(undef, length(gensys2_regimes) - 1)
    Rout  = Vector{Matrix{S}}(undef, length(gensys2_regimes) - 1)
    Cout  = Vector{Vector{S}}(undef, length(gensys2_regimes) - 1)

    # Calculate "uncertain" temporary altpolicy matrices and back out the transition equation
    n_alt = length(prob_vec) - 1
    for (i, reg) in enumerate(gensys2_regimes[2:end - 1])
        # it is assumed T_impl is a vector of the T_{t + 1}^{(temporary altpolicy)} matrices,
        # hence Tbars will be a vector of T_t matrices
        Tbars = prob_vec[1] * (@view systems[1][reg + 1, :TTT][inds, inds]) .+
        sum([prob_vec[j] .* (is_altpol[j - 1] ? (@view systems[j][:TTT][inds, inds]) : (@view systems[j][min(length(systems[j].transitions), reg + 1), :TTT][inds, inds]))
                 for j in 2:n_alt])
        Cbars = prob_vec[1] * (@view systems[1][reg + 1, :CCC][inds]) .+
        sum([prob_vec[j] .* (is_altpol[j - 1] ? (@view systems[j][:CCC][inds]) : (@view systems[j][min(length(systems[j].transitions), reg + 1), :CCC][inds]))
                 for j in 2:n_alt])

        tmp     = Γ2_til * Tbars + Γ0_til
        Tout[i] = tmp \ Γ1_til
        Rout[i] = tmp \ Ψ_til
        Cout[i] = tmp \ (C_til  - Γ2_til * Cbars)
    end

    # Add terminal condition
    l_g2      = last(gensys2_regimes)
    Tout[end] = systems[1][l_g2, :TTT]
    Rout[end] = systems[1][l_g2, :RRR]
    Cout[end] = systems[1][l_g2, :CCC]

    return Tout, Rout, Cout
end

# With time-varying credibility
function gensys2_uncertain_altpol(prob_vec::Vector{Vector{S}}, gensys2_regimes::Vector{Int64}, inds::UnitRange{Int64},
                                  systems::Vector{Union{System, RegimeSwitchingSystem}}, is_altpol::Vector{Bool},
                                  Γ0_til::AbstractMatrix{S}, Γ1_til::AbstractMatrix{S}, Γ2_til::AbstractMatrix{S},
                                  C_til::AbstractVector{S}, Ψ_til::AbstractMatrix{S}) where {S <: Real}


    Tout  = Vector{Matrix{S}}(undef, length(gensys2_regimes) - 1)
    Rout  = Vector{Matrix{S}}(undef, length(gensys2_regimes) - 1)
    Cout  = Vector{Vector{S}}(undef, length(gensys2_regimes) - 1)

    # Calculate "uncertain" temporary altpolicy matrices and back out the transition equation
    n_alt = length(prob_vec[1])
    for (i, reg) in enumerate(gensys2_regimes[2:end - 1]) # gensys2_regimes includes 1 extra regime at beginning for boundary condition
        # it is assumed T_impl is a vector of the T_{t + 1}^{(temporary altpolicy)} matrices,
        # hence Tbars will be a vector of T_t matrices
        Tbars = prob_vec[i][1] * (@view systems[1][reg + 1, :TTT][inds, inds]) .+
        sum([prob_vec[i][j] .* (is_altpol[j - 1] ? (@view systems[j][:TTT][inds, inds]) : (@view systems[j][min(length(systems[j].transitions), reg + 1), :TTT][inds, inds]))
                 for j in 2:n_alt])
        Cbars = prob_vec[i][1] * (@view systems[1][reg + 1, :CCC][inds]) .+
        sum([prob_vec[i][j] .* (is_altpol[j - 1] ? (@view systems[j][:CCC][inds]) : (@view systems[j][min(length(systems[j].transitions), reg + 1), :CCC][inds]))
                 for j in 2:n_alt])

        tmp     = Γ2_til * Tbars + Γ0_til
        Tout[i] = tmp \ Γ1_til
        Rout[i] = tmp \ Ψ_til
        Cout[i] = tmp \ (C_til  - Γ2_til * Cbars)
    end

    # Add terminal condition
    l_g2      = last(gensys2_regimes)
    Tout[end] = systems[1][l_g2, :TTT]
    Rout[end] = systems[1][l_g2, :RRR]
    Cout[end] = systems[1][l_g2, :CCC]

    return Tout, Rout, Cout
end
