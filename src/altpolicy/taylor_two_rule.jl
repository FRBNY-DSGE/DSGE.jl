function taylor_two_rule()
    AltPolicy(:taylor_two_rule, taylor_two_rule_eqcond, taylor_two_rule_solve,
              forecast_init = zero_rate_forecast_init,
              color = RGB(0., 0., 0.5430)) # dark blue
end

"""
```
zero_rate_eqcond(m::AbstractDSGEModel)
```

Solves for the transition equation of `m` under a price level
targeting rule (implemented by adding a price-gap state)
"""
function taylor_two_rule_eqcond(m::AbstractDSGEModel, reg::Int = 1)

    # Get equilibrium condition matrices
    Γ0, Γ1, C, Ψ, Π = if reg <= 11 && reg >= 6
        zero_rate_eqcond(m, reg)
    else
        eqcond(m, reg)
    end
    eq = m.equilibrium_conditions
    endo = m.endogenous_states
    Γ0[eq[:eq_pgap], endo[:pgap_t]] = 1.
    Γ0[eq[:eq_ygap], endo[:ygap_t]] = 1.
    return Γ0, Γ1, C, Ψ, Π
end

"""
```
zero_rate_solve(m::AbstractDSGEModel)
```

Solves for the transition equation of `m` under a price level
targeting rule (implemented by adding a price-gap state)
"""
function taylor_two_rule_solve(m::AbstractDSGEModel; regime_switching::Bool = false, regimes::Vector{Int} = Int[1])
#=    file = JLD2.jldopen("systemtomatch.jld2", "r")
    sys = file["system"]
    n_regs = length(regimes)
    TTTs= Vector{Matrix{Float64}}(undef, n_regs)
    RRRs = Vector{Matrix{Float64}}(undef, n_regs)
    CCCs = Vector{Vector{Float64}}(undef, n_regs)
    for (i, reg) in enumerate(regimes)
        if reg <= 12
            TTTs[i] = sys[reg][:TTT]
            CCCs[i] = sys[reg][:CCC]
            RRRs[i] = sys[reg][:RRR]
        else
            TTTs[i] = sys[12][:TTT]
            CCCs[i] = sys[12][:CCC]
            RRRs[i] = sys[12][:RRR]
        end
    end
    if n_regs == 1
        return TTTs[1], RRRs[1], CCCs[1]
    else
        return TTTs, RRRs, CCCs
    end =#

    m <= Setting(:add_pgap, false)
    m <= Setting(:add_altpolicy_pgap, false)
    m <= Setting(:add_ygap, false)
    m <= Setting(:add_altpolicy_ygap, false)

    # Get equilibrium condition matrices
    if length(regimes) == 1
        Γ0, Γ1, C, Ψ, Π  = taylor_two_rule_eqcond(m, regimes[1])

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
    m <= Setting(:add_pgap, true)
    m <= Setting(:add_altpolicy_pgap, true)
    m <= Setting(:add_ygap, true)
    m <= Setting(:add_altpolicy_ygap, true)
        return TTT, RRR, CCC
    else
        n_regs = get_setting(m, :n_regimes)
        Γ0s = Vector{Matrix{Float64}}(undef, n_regs)
        Γ1s = Vector{Matrix{Float64}}(undef, n_regs)
        Cs = Vector{Vector{Float64}}(undef, n_regs)
        Ψs = Vector{Matrix{Float64}}(undef, n_regs)
        Πs = Vector{Matrix{Float64}}(undef, n_regs)
        for reg in 1:n_regs
            Γ0s[reg], Γ1s[reg], Cs[reg], Ψs[reg], Πs[reg]  = taylor_two_rule_eqcond(m, reg)
        end

        n_regimes = length(regimes)
        TTTs= Vector{Matrix{Float64}}(undef, n_regs)
        RRRs = Vector{Matrix{Float64}}(undef, n_regs)
        CCCs = Vector{Vector{Float64}}(undef, n_regs)

        @show n_regs


        # Solve model
        for reg in 1:n_regs
            if reg <= 5 || reg >= 12
                TTT_gensys, CCC_gensys, RRR_gensys, eu = gensys(Γ0s[reg], Γ1s[reg], Cs[reg], Ψs[reg], Πs[reg], 1+1e-6)

                if !((eu[1] == 1) & (eu[2] == 1))
                    throw(GensysError("Gensys does not give existence"))
                end

                TTT_gensys = real(TTT_gensys)
                RRR_gensys = real(RRR_gensys)
                CCC_gensys = real(CCC_gensys)

                # Augment states
                TTTs[reg], RRRs[reg], CCCs[reg] = DSGE.augment_states(m, TTT_gensys, RRR_gensys, CCC_gensys,reg = reg)
            end
        end
        # solve gensys2 regimes
#        TTT_lift, RRR_lift, CCC_lift = TTT#taylor_rule_solve(m; regime_switching=true, regimes = [12])
        TTT_lift = TTTs[12]
        RRR_lift = RRRs[12]
        CCC_lift = CCCs[12]
        n_endo = length(m.endogenous_states)
        TTT_lift = TTT_lift[1:n_endo, 1:n_endo]
        RRR_lift = RRR_lift[1:n_endo, :]
        CCC_lift = CCC_lift[1:n_endo]

        gensys2_regimes = collect(5:12)
        Tcal, Rcal, Ccal = gensys2(m, Γ0s[gensys2_regimes], Γ1s[gensys2_regimes], Cs[gensys2_regimes], Ψs[gensys2_regimes], Πs[gensys2_regimes],
                                   TTT_lift, RRR_lift, CCC_lift,
                                   T_switch = length(gensys2_regimes)-1)
        Tcal[end] = TTT_lift
        Rcal[end] = RRR_lift
        Ccal[end] = CCC_lift

        populate_reg = gensys2_regimes[2:end]
        for (i, reg) in enumerate(populate_reg)
            TTTs[reg], RRRs[reg], CCCs[reg] = augment_states(m, Tcal[i], Rcal[i], Ccal[i], reg=reg)
        end
    end
    JLD2.jldopen("hybridsystem.jld2", true, true, true, IOStream) do file
        write(file, "TTTs", TTTs)
        write(file, "CCCs", CCCs)
    end

    m <= Setting(:add_pgap, true)
    m <= Setting(:add_altpolicy_pgap, true)
    m <= Setting(:add_ygap, true)
    m <= Setting(:add_altpolicy_ygap, true)

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
