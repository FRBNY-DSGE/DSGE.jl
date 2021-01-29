function zero_rate()
    AltPolicy(:zero_rate, zero_rate_eqcond, zero_rate_solve,
              forecast_init = zero_rate_forecast_init,
              color = RGB(0., 0., 0.5430)) # dark blue
end


function zero_rate_replace_eq_entries(m::AbstractDSGEModel,
                                      Γ0::Matrix{Float64}, Γ1::Matrix{Float64},
                                      C::Vector{Float64}, Ψ::Matrix{Float64}, Π::Matrix{Float64})

    eq             = m.equilibrium_conditions
    endo           = m.endogenous_states

    Γ0[eq[:eq_mp], :] .= 0.
    Γ1[eq[:eq_mp], :] .= 0.
    Γ0[eq[:eq_mp], endo[:R_t]] = 1.0
    C[eq[:eq_mp]] = 0.1 / 4. - m[:Rstarn]

    return Γ0, Γ1, C, Ψ, Π
end

"""
```
zero_rate_eqcond(m::AbstractDSGEModel)
```

Solves for the transition equation of `m` under a price level
targeting rule (implemented by adding a price-gap state)
"""
function zero_rate_eqcond(m::AbstractDSGEModel, reg::Int = 1)

    # Get equilibrium condition matrices
    Γ0, Γ1, C, Ψ, Π = eqcond(m, reg)
    Γ0, Γ1, C, Ψ, Π = zero_rate_replace_eq_entries(m, Γ0, Γ1, C, Ψ, Π)

    return Γ0, Γ1, C, Ψ, Π
end

"""
```
zero_rate_solve(m::AbstractDSGEModel)
```

Solves for the transition equation of `m` under a price level
targeting rule (implemented by adding a price-gap state)
"""
function zero_rate_solve(m::AbstractDSGEModel; regime_switching::Bool = false, regimes::Vector{Int} = Int[1])

    # Get equilibrium condition matrices
    if length(regimes) == 1
        Γ0, Γ1, C, Ψ, Π  = zero_rate_eqcond(m, regimes[1])

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
            Γ0s[reg], Γ1s[reg], Cs[reg], Ψs[reg], Πs[reg]  = zero_rate_eqcond(m, reg)
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
    end

    return TTTs, RRRs, CCCs
end

"""
```
init_zero_rate_forecast(m::AbstractDSGEModel, shocks::Matrix{T}, final_state::Vector{T})
```

Adjust shocks matrix and final state vector for forecasting under the ZERO_RATE rule
"""
function zero_rate_forecast_init(m::AbstractDSGEModel, shocks::Matrix{T}, final_state::Vector{T};
                                 cond_type::Symbol = :none) where {T<:AbstractFloat}
    return shocks, final_state
end
