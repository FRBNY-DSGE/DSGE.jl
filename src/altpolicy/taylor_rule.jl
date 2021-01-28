"""
```
taylor_rule()
```

creates a Taylor rule with monetary policy shocks and coefficients on the
contemporaneous inflation gap, output gap, and growth in the output gap.
"""
function taylor_rule()
    AltPolicy(:taylor_rule, taylor_rule_eqcond, taylor_rule_solve,
              forecast_init = taylor_rule_forecast_init,
              color = RGB(0., 0., 0.5430)) # dark blue
end

"""
```
taylor_replace_eq_entries(m::AbstractDSGEModel,
                          Γ0::Matrix{Float64}, Γ1::Matrix{Float64},
                          C::Vector{Float64}, Ψ::Matrix{Float64}, Π::Matrix{Float64})taylor_replace_eq_entries()
```

implements a Taylor rule with monetary policy shocks and coefficients on the
contemporaneous inflation gap, output gap, and growth in the output gap.
"""
function taylor_rule_replace_eq_entries(m::AbstractDSGEModel,
                                        Γ0::Matrix{Float64}, Γ1::Matrix{Float64},
                                        C::Vector{Float64}, Ψ::Matrix{Float64}, Π::Matrix{Float64})

    # add law of motion for pgap
    eq        = m.equilibrium_conditions
    endo      = m.endogenous_states

    # Zero out old policy rule
    Γ0[eq[:eq_mp], :]              .= 0.
    Γ1[eq[:eq_mp], :]              .= 0.
    Ψ[eq[:eq_mp], :]               .= 0.

    # Replace monetary policy rule
    Γ0[eq[:eq_mp], endo[:R_t]]      = 1.
    Γ1[eq[:eq_mp], endo[:R_t]]      = m[:ρ]
    C[ eq[:eq_mp]]                  = 0.


    Γ0[eq[:eq_mp], endo[:π_t]]      = -(1. - m[:ρ]) * m[:ψ1]
    Γ0[eq[:eq_mp], endo[:π_star_t]] = (1. - m[:ρ]) * m[:ψ1]
    Γ0[eq[:eq_mp], endo[:y_t]]      = -(1. - m[:ρ]) * m[:ψ2] - m[:ψ3]
    Γ0[eq[:eq_mp], endo[:y_f_t]]    = (1. - m[:ρ]) * m[:ψ2] + m[:ψ3]
    Γ1[eq[:eq_mp], endo[:y_t]]      = -m[:ψ3]
    Γ1[eq[:eq_mp], endo[:y_f_t]]    = m[:ψ3]

    # Add MP shocks
    Γ0[eq[:eq_mp], endo[:rm_t]]     = -1.

    return Γ0, Γ1, C, Ψ, Π
end

"""
```
taylor_rule_eqcond(m::AbstractDSGEModel)
```

Solves for the transition equation of `m` under a Taylor Rule.
"""
function taylor_rule_eqcond(m::AbstractDSGEModel, reg::Int = 1)

    # Get equilibrium condition matrices
    Γ0, Γ1, C, Ψ, Π = eqcond(m, reg)

    # Update parameters to regime `reg` since `eqcond` will reset to the first regime
    # Note that this step is necessary (in general) since
    # taylor_replace_eq_entries will use parameters as coefficients for
    # the monetary policy rule
    for para in m.parameters
        if !isempty(para.regimes)
            if (haskey(get_settings(m), :model2para_regimes) ? haskey(get_setting(m, :model2para_regime), para.key) : false)
                ModelConstructors.toggle_regime!(para, reg, get_setting(m, :model2para_regime)[para.key])
            else
                ModelConstructors.toggle_regime!(para, reg)
            end
        end
    end

    Γ0, Γ1, C, Ψ, Π = taylor_rule_replace_eq_entries(m, Γ0, Γ1, C, Ψ, Π)

    # Reset to the first regime
    ModelConstructors.toggle_regime!(m.parameters, 1)

    return Γ0, Γ1, C, Ψ, Π
end

"""
```
taylor_rule_solve(m::AbstractDSGEModel)
```

Solves for the transition equation of `m` under a Taylor Rule.
"""
function taylor_rule_solve(m::AbstractDSGEModel; regime_switching::Bool = false, regimes::Vector{Int} = Int[1])

    # Get equilibrium condition matrices
    if length(regimes) == 1
        Γ0, Γ1, C, Ψ, Π  = taylor_rule_eqcond(m, regimes[1])
        TTT_gensys, CCC_gensys, RRR_gensys, eu = gensys(Γ0, Γ1, C, Ψ, Π, 1+1e-6, verbose = :low)

        # Check for LAPACK exception, existence and uniqueness
        if eu[1] != 1 || eu[2] != 1
            throw(GensysError())
        end

        TTT_gensys = real(TTT_gensys)
        RRR_gensys = real(RRR_gensys)
        CCC_gensys = real(CCC_gensys)

        # Augment states
        TTT, RRR, CCC = DSGE.augment_states(m, TTT_gensys, RRR_gensys, CCC_gensys; regime_switching = regime_switching,
                                            reg = regimes[1])
        return TTT, RRR, CCC
    else
        Γ0s = Vector{Matrix{Float64}}(undef, length(regimes))
        Γ1s = Vector{Matrix{Float64}}(undef, length(regimes))
        Cs = Vector{Vector{Float64}}(undef, length(regimes))
        Ψs = Vector{Matrix{Float64}}(undef, length(regimes))
        Πs = Vector{Matrix{Float64}}(undef, length(regimes))
        for reg in regimes
            Γ0s[reg], Γ1s[reg], Cs[reg], Ψs[reg], Πs[reg]  = taylor_rule_eqcond(m, reg)
        end

        n_regimes = length(regimes)
        TTTs= Vector{Matrix{Float64}}(undef, n_regimes)
        RRRs = Vector{Matrix{Float64}}(undef, n_regimes)
        CCCs = Vector{Vector{Float64}}(undef, n_regimes)

        # Solve model
        for reg in regimes
            TTT_gensys, CCC_gensys, RRR_gensys, eu = gensys(Γ0s[reg], Γ1s[reg], Cs[reg], Ψs[reg], Πs[reg], 1+1e-6)

            if !((eu[1] == 1) & (eu[2] == 1))
                throw(GensysError("Gensys does not give existence"))
            end
            TTT_gensys = real(TTT_gensys)
            RRR_gensys = real(RRR_gensys)
            CCC_gensys = reshape(CCC_gensys, size(CCC_gensys, 1))

            # Augment states
            TTTs[reg], RRRs[reg], CCCs[reg] = DSGE.augment_states(m, TTT_gensys, RRR_gensys, CCC_gensys;
                                                             regime_switching = regime_switching,
                                                             reg = reg)
        end
        return TTTs, RRRs, CCCs
    end
end

"""
```
init_taylor_rule_forecast(m::AbstractDSGEModel, shocks::Matrix{T}, final_state::Vector{T})
```

Adjust shocks matrix and final state vector for forecasting under a Taylor Rule.
"""
function taylor_rule_forecast_init(m::AbstractDSGEModel, shocks::Matrix{T}, final_state::Vector{T};
                                   cond_type::Symbol = :none) where {T<:AbstractFloat}
    return shocks, final_state
end
