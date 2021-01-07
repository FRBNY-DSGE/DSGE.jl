function zero_rate()
    AltPolicy(:zero_rate, zero_rate_eqcond, zero_rate_solve,
              replace_eqcond = zero_rate_replace_eq_entries,
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
    C[eq[:eq_mp]] = 0.0 / 4. - m[:Rstarn]

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

    # get the old indices
    old_states = sort!(collect(values(m.endogenous_states)))
    old_eqs    = sort!(collect(values(m.equilibrium_conditions)))

    # Get equilibrium condition matrices
    Γ0, Γ1, C, Ψ, Π  = eqcond(m, reg)

    for para in m.parameters
        if !isempty(para.regimes)
            if (haskey(get_settings(m), :model2para_regime) ? haskey(get_setting(m, :model2para_regime), para.key) : false)
                ModelConstructors.toggle_regime!(para, reg, get_setting(m, :model2para_regime)[para.key])
            else
                ModelConstructors.toggle_regime!(para, reg)
            end
        end
    end

    Γ0, Γ1, C, Ψ, Π = zero_rate_replace_eq_entries(m, Γ0, Γ1, C, Ψ, Π)

    for para in m.parameters
        if !isempty(para.regimes)
            ModelConstructors.toggle_regime!(para, 1)
        end
    end

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
    return DSGE.altpolicy_solve(m, zero_rate(); regime_switching = regime_switching, regimes = regimes)
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
