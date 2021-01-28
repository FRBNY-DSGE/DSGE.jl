"""
```
ait()
```

creates a permanent alternative policy of "Average Inflation Targeting".
To run this permanent policy, a user needs to run the following code:

```
setup_permanent_altpolicy!(m, ait())
m <= Setting(:pgap_type, :ait)
m <= Setting(:pgap_value, 0.)
```

The first line initializes the alternative policy, and the next two lines
indicates what information is represented by the price gap (`pgap`). It is
standard to set `pgap_value = 0.`, which means that the AIT rule
implements an initial price gap target of 0 percent.

A user should also check that temporary alternative policies are turned off,
namely that the settings `replace_eqcond` and `gensys2` should be set to `false`.
"""
function ait()
    AltPolicy(:ait, ait_eqcond, ait_solve,
              forecast_init = ait_forecast_init,
              color = RGB(0., 0., 0.5430)) # dark blue
end

function ait_replace_eq_entries(m::AbstractDSGEModel,
                                Γ0::Matrix{Float64}, Γ1::Matrix{Float64},
                                C::Vector{Float64}, Ψ::Matrix{Float64}, Π::Matrix{Float64})

    # add law of motion for pgap
    eq             = m.equilibrium_conditions
    endo           = m.endogenous_states

    Thalf = haskey(get_settings(m), :ait_Thalf) ? get_setting(m, :ait_Thalf) : 10.
    ρ_ait = exp(log(0.5) / Thalf)
    ρ     = 0.0
    φ     = 0.25

    # This assumes that the inflation target is the model's steady state
    Γ0[eq[:eq_pgap], endo[:pgap_t]]  =  1.
    Γ0[eq[:eq_pgap], endo[:π_t]]     = -1.
    Γ1[eq[:eq_pgap], endo[:pgap_t]]  = ρ_ait

###
#= Γ0[eq[:eq_pgap], endo[:y_t]]     = -1.
    Γ0[eq[:eq_pgap], endo[:z_t]]     = -1.
    Γ1[eq[:eq_pgap], endo[:y_t]]     =  -1.=#
###

    # Zero out old policy rule
    Γ0[eq[:eq_mp], :]             .= 0.
    Γ1[eq[:eq_mp], :]             .= 0.
    Ψ[eq[:eq_mp], :]              .= 0.

    # replace monetary policy rule
    Γ0[eq[:eq_mp], endo[:R_t]]       = 1.

    Γ0[eq[:eq_mp], endo[:pgap_t]]    = -φ*(1/(1-ρ_ait))
    Γ1[eq[:eq_mp], endo[:R_t]]       = ρ
    C[eq[:eq_mp]]                    = 0.

    return Γ0, Γ1, C, Ψ, Π
end

"""
```
ait_eqcond(m::AbstractDSGEModel)
```

Solves for the transition equation of `m` under a price level
targeting rule (implemented by adding a price-gap state)
"""
function ait_eqcond(m::AbstractDSGEModel, reg::Int = 1)

    # get the old indices
    old_states = sort!(collect(values(m.endogenous_states)))
    old_eqs    = sort!(collect(values(m.equilibrium_conditions)))

    # Get equilibrium condition matrices
    Γ0_noaltpol, Γ1_noaltpol, C_noaltpol, Ψ_noaltpol, Π_noaltpol  = eqcond(m, reg)

    # check to make sure this state hasn't already been added, then:
    # 1. add the pgap state to the model
    # 2. increment the indices of every augmented state

    if !in(:pgap_t, keys(m.endogenous_states))
        # add pgap state to endogenous_staets
        m.endogenous_states[:pgap_t]       = n_states(m) + 1
        m.equilibrium_conditions[:eq_pgap] = n_equilibrium_conditions(m) + 1
        # now n_states(m) and n_equilibrium_conditions(m) will calculate appropriately

        # increment every state in endogenous_states_augmented
        aug = m.endogenous_states_augmented
        map(x -> aug[x] = aug[x]+1, collect(keys(aug)))
    end

    nstates = n_states(m)

    # fill in new Γ0, Γ1, C, Ψ, Π
    Γ0 = zeros(Float64, nstates, nstates)
    Γ0[old_eqs, old_states] = Γ0_noaltpol

    Γ1 = zeros(Float64, nstates, nstates)
    Γ1[old_eqs, old_states] = Γ1_noaltpol

    C  = zeros(Float64, nstates)
    C[old_eqs, :] = C_noaltpol

    Ψ  = zeros(Float64, nstates, n_shocks_exogenous(m))
    Ψ[old_eqs, :] = Ψ_noaltpol

    Π  = zeros(Float64, nstates, n_shocks_expectational(m))
    Π[old_eqs, :] = Π_noaltpol

    for para in m.parameters
        if !isempty(para.regimes)
            if (haskey(get_settings(m), :model2para_regime) ? haskey(get_setting(m, :model2para_regime), para.key) : false)
                ModelConstructors.toggle_regime!(para, reg, get_setting(m, :model2para_regime)[para.key])
            else
                ModelConstructors.toggle_regime!(para, reg)
            end
        end
    end

    Γ0, Γ1, C, Ψ, Π = ait_replace_eq_entries(m, Γ0, Γ1, C, Ψ, Π)
    for para in m.parameters
        if !isempty(para.regimes)
            ModelConstructors.toggle_regime!(para, 1)
        end
    end

    return Γ0, Γ1, C, Ψ, Π
end

"""
```
ait_solve(m::AbstractDSGEModel)
```

Solves for the transition equation of `m` under a price level
targeting rule (implemented by adding a price-gap state, which uses the value
held in the setting `get_setting(m, :pgap_value)`. A value of 0. is standard.)
"""
function ait_solve(m::AbstractDSGEModel; regime_switching::Bool = false, regimes::Vector{Int} = Int[1])

    # Get equilibrium condition matrices
    if length(regimes) == 1
        Γ0, Γ1, C, Ψ, Π  = ait_eqcond(m, regimes[1])
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
            Γ0s[reg], Γ1s[reg], Cs[reg], Ψs[reg], Πs[reg]  = ait_eqcond(m, reg)
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
init_ait_forecast(m::AbstractDSGEModel, shocks::Matrix{T}, final_state::Vector{T})
```

Adjust shocks matrix and final state vector for forecasting under the AIT rule
"""
function ait_forecast_init(m::AbstractDSGEModel, shocks::Matrix{T}, final_state::Vector{T};
                           cond_type::Symbol = :none) where {T<:AbstractFloat}

    pgap_t = m.endogenous_states[:pgap_t]
    # THIS IS WRONG --> #final_state[pgap_t] = -3.5 # todo: figure out how to program automatically
   final_state = vcat(final_state[1:pgap_t-1], [-get_setting(m, :pgap_value)], final_state[pgap_t+1:end]) #-3.5

    return shocks, final_state
end
