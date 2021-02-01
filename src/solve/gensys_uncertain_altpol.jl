function gensys_uncertain_altpol(m::AbstractDSGEModel, prob_vec::AbstractVector{S},
                                 altpolicies::Vector{AltPolicy} = [default_policy()];
                                 regime_switching::Bool = false, regimes::Vector{Int} = Int[1],
                                 TTT::Matrix{S} = Matrix{S}(undef, 0, 0),
                                 Γ0s::Vector{Matrix{S}} = Vector{Matrix{S}}(undef, 0),
                                 Γ1s::Vector{Matrix{S}} = Vector{Matrix{S}}(undef, 0),
                                 Cs::Vector{Vector{S}} = Vector{Vector{S}}(undef, 0),
                                 Ψs::Vector{Matrix{S}} = Vector{Matrix{S}}(undef, 0),
                                 Πs::Vector{Matrix{S}} = Vector{Matrix{S}}(undef, 0)) where {S <: Real}

    @assert sum(prob_vec) == 1. "The vector of probabilities must sum to 1"

    @assert length(altpolicies) == length(prob_vec) - 1

    altpolicies = altpolicies[prob_vec[2:end] .> 0.0]
    prob_vec = vcat(prob_vec[1], prob_vec[2:end][prob_vec[2:end] .> 0.0])
    if regime_switching
        if length(regimes) == 1
            Γ0, Γ1, C, Ψ, Π = haskey(get_settings(m), :regime_eqcond_info) &&
                haskey(get_setting(m, :regime_eqcond_info), regimes[1]) ?
                get_setting(m, :regime_eqcond_info)[regimes[1]].alternative_policy.eqcond(m, regimes[1]) :
                eqcond(m, regimes[1])

            if isempty(TTT)
                assert_cond = haskey(get_settings(m), :uncertain_altpolicy) ? !get_setting(m, :uncertain_altpolicy) : true
                @assert assert_cond "The Setting :uncertain_altpolicy must be false if TTT is empty."
                TTT, _, _ = solve(m; regime_switching = regime_switching, regimes = regimes)
            end

            Γ0_til, Γ1_til, Γ2_til, C_til, Ψ_til = gensys_to_predictable_form(Γ0, Γ1, C, Ψ, Π)

            inds = 1:n_states(m) # Don't want augmented states
            Td = Vector{Matrix{S}}(undef, length(prob_vec))
            Cd = Vector{Vector{S}}(undef, length(prob_vec))
            Td[1] = TTT[inds, inds]
            Cd[1] = C[inds]

            # Turn off uncertain altpolicy b/c when computing the state space system under alternative policies,
            # we want to treat these alternative policies as perfectly credible.
            is_uncertain_altpol = haskey(get_settings(m), :uncertain_altpolicy) && get_setting(m, :uncertain_altpolicy)
            if is_uncertain_altpol
                m <= Setting(:uncertain_altpolicy, false)
            end

            # Solve each alternative policy
            for (i, altpolicy) in enumerate(altpolicies)
                j = i + 1
                tmpT, _, tmpC = if altpolicy.solve == solve      # if solve is being called for gensys_uncertain_altpol
                    altpolicy.solve(m; regime_switching = false) # then it is interpreted as the "default" policy, and
                else                                             # it is further assumed that we can safely ignore regime-switching
                    altpolicy.solve(m; regime_switching = regime_switching, regimes = regimes)
                end
                Td[j] = tmpT[inds, inds]
                Cd[j] = tmpC[inds]
            end

            # Turn uncertain_altpolicy back on
            if is_uncertain_altpol
                m <= Setting(:uncertain_altpolicy, true)
            end

            T̄ = sum([p .* Td[i] for (i, p) in enumerate(prob_vec)])
            C̄ = sum([p .* Cd[i] for (i, p) in enumerate(prob_vec)])
            Lmat = (Γ2_til * T̄ + Γ0_til)
            Tcal = Lmat \ Γ1_til
            Rcal = Lmat \ Ψ_til
            Ccal = Lmat \ (C_til - Γ2_til * C̄)
            return Tcal, Rcal, Ccal
        else
            error("Regime switching for multiple regimes has not been completed for gensys_uncertain_altpol")
            if isempty(Γ0s) || isempty(Γ1s) || isempty(Cs) || isempty(Ψs) || isempty(Πs)
                # If any of these are empty, we recompute the relevant matrices
                for fcast_reg in regimes
                    Γ0s[fcast_reg], Γ1s[fcast_reg], Cs[fcast_reg], Ψs[fcast_reg], Πs[fcast_reg] =
                        haskey(get_settings(m), :regime_eqcond_info) ?
                        (haskey(get_setting(m, :regime_eqcond_info), fcast_reg) ?
                         get_setting(m, :regime_eqcond_info)[fcast_reg].alternative_policy.eqcond(m, fcast_reg) :
                         eqcond(m, fcast_reg)) : eqcond(m, fcast_reg)
                end
            end

            # Solve model for creating the "default" matrices in each period
            TTTs_gensys, RRRs_gensys, CCCs_gensys, eus =
                map((Γ0, Γ1, C, Ψ, Π) -> gensys(Γ0, Γ1, C, Ψ, Π, 1 + 1e-6, verbose = verbose), Γ0s, Γ1s, Cs, Ψs, Πs)

            # Check for LAPACK exception, existence and uniqueness
            for eu in eus
                if eu[1] != 1 || eu[2] != 1
                    throw(GensysError())
                end
            end

            TTTs_gensys = map(TTT_gensys -> real(TTT_gensys), TTTs_gensys)
            RRRs_gensys = map(RRR_gensys -> real(RRR_gensys), RRRs_gensys)
            CCCs_gensys = map(CCC_gensys -> real(CCC_gensys), CCCs_gensys)

            Γ0_tils, Γ1_tils, Γ2_tils, C_tils, Ψ_tils =
                map((Γ0, Γ1, C, Ψ, Π) -> gensys_to_predictable_form(Γ0, Γ1, C, Ψ, Π), Γ0s, Γ1s, Cs, Ψs, Πs)

            # Rows are different altpolicies in each regime, columns are different regimes
            Td = Matrix{Matrix{S}}(undef, length(prob_vec))
            Cd = Matrix{Vector{S}}(undef, length(prob_vec))
            Td[1, :] = TTTs_gensys
            Cd[1, :] = CCCs_gensys
            for fcast_reg in regimes
                for (i, altpolicy) in enumerate(altpolicies)
                    j = i + 1
                    tmpT, _, tmpC = altpolicy.solve(m; regime_switching = regime_switching, regimes = regimes)
                    Td[j, fcast_reg] = tmpT
                    Cd[j, fcast_reg] = tmpC
                end
            end

            T̄s = map(reg -> sum([p .* Td_reg[j, reg] for (j, p) in enumerate(prob_vec)]), regimes)
            C̄s = map(reg -> sum([p .* Cd_reg[j, reg] for (j, p) in enumerate(prob_vec)]), regimes)

            Tcal = Vector{Matrix{S}}(undef, length(regimes))
            Rcal = Vector{Matrix{S}}(undef, length(regimes))
            Ccal = Vector{Vector{S}}(undef, length(regimes))

            for i in length(regimes)
                Lmat = (Γ2_tils[i] * T̄s[i] + Γ0_tils[i])
                Tcals[i] = Lmat \ Γ1_tils[i]
                Rcals[i] = Lmat \ Ψ_tils[i]
                Ccals[i] = Lmat \ (C_tils[i] - Γ2_tils[i] * C̄s[i])
            end

            return Tcals, Rcals, Ccals
        end
    else
        Γ0, Γ1, C, Ψ, Π = haskey(get_settings(m), :regime_eqcond_info) ?
            (haskey(get_setting(m, :regime_eqcond_info), regimes[1]) ?
             get_setting(m, :regime_eqcond_info)[regimes[1]].alternative_policy.eqcond(m, regimes[1]) :
             eqcond(m, regimes[1])) : eqcond(m, regimes[1])

        if isempty(TTT)
            assert_cond = haskey(get_settings(m), :uncertain_altpolicy) ? !get_setting(m, :uncertain_altpolicy) : true
            @assert assert_cond "The Setting :uncertain_altpolicy must be false if TTT is empty."
            TTT, _, _ = solve(m)
        end

        Γ0_til, Γ1_til, Γ2_til, C_til, Ψ_til = gensys_to_predictable_form(Γ0, Γ1, C, Ψ, Π)

        # Turn off uncertain altpolicy b/c when computing the state space system under alternative policies,
        # we want to treat these alternative policies as perfectly credible.
        is_uncertain_altpol = haskey(get_settings(m), :uncertain_altpolicy) && get_setting(m, :uncertain_altpolicy)
        if is_uncertain_altpol
            m <= Setting(:uncertain_altpolicy, false)
        end

        inds = 1:n_states(m) # Don't want augmented states
        Td = Vector{Matrix{S}}(undef, length(prob_vec))
        Cd = Vector{Vector{S}}(undef, length(prob_vec))
        Td[1] = TTT[inds, inds]
        Cd[1] = C[inds]
        for (i, altpolicy) in enumerate(altpolicies)
            j = i + 1
            tmpT, _, tmpC = altpolicy.solve(m) # regime_switching is false, so it is assumed that it is unnecessary here
            Td[j] = tmpT[inds, inds]
            Cd[j] = tmpC[inds]
        end

        # Turn uncertain_altpolicy back on
        if is_uncertain_altpol
            m <= Setting(:uncertain_altpolicy, true)
        end

        T̄ = sum([p .* Td[i] for (i, p) in enumerate(prob_vec)])
        C̄ = sum([p .* Cd[i] for (i, p) in enumerate(prob_vec)])

        Lmat = (Γ2_til * T̄ + Γ0_til)
        Tcal = Lmat \ Γ1_til
        Rcal = Lmat \ Ψ_til
        Ccal = Lmat \ (C_til - Γ2_til * C̄)

        return Tcal, Rcal, Ccal
    end
end

function gensys_uncertain_altpol(m::AbstractDSGEModel, prob_vec::AbstractVector{S}, regime::Int, TTT::Matrix{S},
                                 TTTs_alt::Vector{Matrix{S}}, CCCs_alt::Vector{Vector{S}},
                                 is_altpol::Vector{Bool}) where {S <: Real}

    @assert sum(prob_vec) == 1. "The vector of probabilities must sum to 1"

    @assert length(TTTs_alt) == length(prob_vec) - 1

    has_pos_prob = findall((@view prob_vec[2:end]) .> 0.) # Note that has_pos_prob has 1 less length than prob_vec

    Γ0, Γ1, C, Ψ, Π = haskey(get_setting(m, :regime_eqcond_info), regime) ?
        get_setting(m, :regime_eqcond_info)[regime].alternative_policy.eqcond(m, regime) :
            eqcond(m, regime)

    Γ0_til, Γ1_til, Γ2_til, C_til, Ψ_til = gensys_to_predictable_form(Γ0, Γ1, C, Ψ, Π)

    T̅, C̅ = prob_vec[1] == 0. ? (zero(TTT), zero(C)) : (prob_vec[1] .* TTT, prob_vec[1] .* C)
    for i in has_pos_prob
        if is_altpol[i]
            T̅ .+= prob_vec[i + 1] * TTTs_alt[i]
            C̅ .+= prob_vec[i + 1] * CCCs_alt[i]
        else
            T̅ .+= prob_vec[i + 1] * TTTs_alt[i][regime]
            C̅ .+= prob_vec[i + 1] * CCCs_alt[i][regime]
        end
    end

    Lmat = (Γ2_til * T̅ + Γ0_til)
    Tcal = Lmat \ Γ1_til
    Rcal = Lmat \ Ψ_til
    Ccal = Lmat \ (C_til - Γ2_til * C̅)

    return Tcal, Rcal, Ccal
end
