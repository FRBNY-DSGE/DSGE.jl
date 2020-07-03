function gensys_prob_regime_switch(m::AbstractDSGEModel, prob_vec::AbstractVector{S},
                                   altpolicies::Vector{AltPolicy} = [get_setting(m, :alternative_policy)];
                                   apply_altpolicy::Bool = false, regime_switching::Bool = false,
                                   regimes::Union{Int, Vector{Int}, UnitRange{Int}} = 1,
                                   Γ0s::Vector{Matrix{S}} = Vector{Matrix{S}}(undef, 0),
                                   Γ1s::Vector{Matrix{S}} = Vector{Matrix{S}}(undef, 0),
                                   Cs::Vector{Vector{S}} = Vector{Vector{S}}(undef, 0),
                                   Ψs::Vector{Matrix{S}} = Vector{Matrix{S}}(undef, 0),
                                   Πs::Vector{Matrix{S}} = Vector{Matrix{S}}(undef, 0)) where {S <: Real}

    @assert sum(prob_vec) == 1. "The vector of probabilities must sum to 1"
    # prob[1] = probability of usual rule

    @assert length(altpolicies) == length(prob_vec) - 1

    if regime_switching
        if isempty(Γ0s) || isempty(Γ1s) || isempty(Cs) || isempty(Ψs) || isempty(Πs)
            # If any of these are empty, we recompute the relevant matrices
            for fcast_reg in regimes
                Γ0s[fcast_reg], Γ1s[fcast_reg], Cs[fcast_reg], Ψs[fcast_reg], Πs[fcast_reg] = eqcond(m, fcast_reg)
            end

            # ALSO NEED TO ADD THE APPLY ALTPOLICY OPTION
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
                j = i + 1 # may refactor to directly use solve(m; apply_altpolicy = true) by updating the alternative policy
                tmpT, _, tmpC = altpolicy.solve(m) # no regime switching for now ; regime_switching = regime_switching, regimes = regimes)
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
    else
        Γ0, Γ1, C, Ψ, Π = if apply_altpolicy
            get_setting(m, :alternative_policy).eqcond(m)
        else
            eqcond(m)
        end

        T, R, _ = solve(m; apply_altpolicy = apply_altpolicy)
        # May need to update later to have m's policy also have apply_altpolicy = true, with the understanding that
        # in the standard case, the altpolicy should be the "historical" altpolicy, i.e. not really altpolicy
        Γ0_til, Γ1_til, Γ2_til, C_til, Ψ_til = gensys_to_predictable_form(Γ0, Γ1, C, Ψ, Π)

        inds = 1:n_states(m) # Don't want augmented states
        Td = Vector{Matrix{S}}(undef, length(prob_vec))
        Cd = Vector{Vector{S}}(undef, length(prob_vec))
        Td[1] = T[inds, inds]
        Cd[1] = C[inds]
        for (i, altpolicy) in enumerate(altpolicies)
            j = i + 1 # may refactor to directly use solve(m; apply_altpolicy = true) by updating the alternative policy
            tmpT, _, tmpC = altpolicy.solve(m) # no regime switching for now ; regime_switching = regime_switching, regimes = regimes)
            Td[j] = tmpT[inds, inds]
            Cd[j] = tmpC[inds]
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

#=
# Uncomment these lines to test
using DSGE, ModelConstructors, Dates, OrderedCollections, Test
m = Model1002("ss10"; custom_settings = Dict{Symbol, Setting}(:add_pgap => Setting(:add_pgap, true)))
m <= Settng(:pgap_value, 12.)
histpolicy = get_setting(m, :alternative_policy)
altpolicies = DSGE.AltPolicy[DSGE.ngdp()]
probs = Float64[.05, .95]

T, R, C = solve(m)
sys = compute_system(m) # get the measurement matrices
m <= Setting(:alternative_policy, DSGE.ngdp())
Tngdp, Rngdp, Cngdp = solve(m; apply_altpolicy = true)
m <= Setting(:alternative_policy, histpolicy)

# Test with fake weights
T1, R1, C1 = gensys_prob_regime_switch(m, [1., 0.], altpolicies)
inds = 1:n_states(m)
@test T1 ≈ T[inds, inds]
@test R1 ≈ R[inds, :]
@test C1 ≈ C[inds]

@test_throws AssertionError gensys_prob_regime_switch(m, [1., .0001], altpolicies)

m <= Setting(:alternative_policy, DSGE.ngdp())
Tngdp1, Rngdp1, Cngdp1 = gensys_prob_regime_switch(m, [1., 0.], DSGE.AltPolicy[histpolicy]; apply_altpolicy = true)
@test Tngdp1 ≈ Tngdp[inds, inds]
@test Rngdp1 ≈ Rngdp[inds, :]
@test Cngdp1 ≈ Cngdp[inds]

# Test 2
m <= Setting(:alternative_policy, histpolicy)
Tcal, Rcal, Ccal = gensys_prob_regime_switch(m, [0.95, 0.05], altpolicies)
Tcal, Rcal, Ccal = DSGE.augment_states(m, Tcal, Rcal, Ccal)
# Tusual, Rusual, Cusual = altpolicies[1].solve(m)
Tusual, Rusual, Cusual = solve(m)
usual_sys = deepcopy(sys)
cal_sys  = deepcopy(sys)
usual_sys.transition.TTT[:, :] = Tusual
usual_sys.transition.RRR[:, :] = Rusual
usual_sys.transition.CCC[:, :] = Cusual
cal_sys.transition.TTT[:, :] = Tcal
cal_sys.transition.RRR[:, :] = Rcal
cal_sys.transition.CCC[:, :] = Ccal
init_state = zeros(size(Tcal, 1))
h = 40
# Random.seed!(2020)
# μ = zeros(n_shocks_exogenous(m))
# σ = sqrt.(cal_sys[:QQ])
# dist = DegenerateMvNormal(μ, σ)
# shocks = rand(dist, h)
# shocks[:, 2:end] .= 0.
shocks = zeros(n_shocks_exogenous(m), h)
shocks[m.exogenous_shocks[:b_sh], 1] = -.05


_, usual_obs, _, _ = forecast(usual_sys, init_state, shocks)
_, cal_obs, _, _ = forecast(cal_sys, init_state, shocks)

using Plots
p2 = plot(1:h, usual_obs[m.observables[:obs_gdp], :], label = "Usual Policy", linewidth = 3)
plot!(1:h, cal_obs[m.observables[:obs_gdp], :], label = "Regime-Switching Usual, P[Usual] = 0.95",
      legend = :bottomright, linewidth = 3)

# Test 3
m <= Setting(:alternative_policy, DSGE.ngdp())
Tcal, Rcal, Ccal = gensys_prob_regime_switch(m, [0.95, 0.05], AltPolicy[histpolicy]; apply_altpolicy = true)
Tcal, Rcal, Ccal = DSGE.augment_states(m, Tcal, Rcal, Ccal)
Tngdp, Rngdp, Cngdp = DSGE.ngdp().solve(m)
ngdp_sys = deepcopy(sys)
cal_sys  = deepcopy(sys)
ngdp_sys.transition.TTT[:, :] = Tngdp
ngdp_sys.transition.RRR[:, :] = Rngdp
ngdp_sys.transition.CCC[:, :] = Cngdp
cal_sys.transition.TTT[:, :] = Tcal
cal_sys.transition.RRR[:, :] = Rcal
cal_sys.transition.CCC[:, :] = Ccal
init_state = zeros(size(Tcal, 1))
h = 40
# Random.seed!(2020)
# μ = zeros(n_shocks_exogenous(m))
# σ = sqrt.(cal_sys[:QQ])
# dist = DegenerateMvNormal(μ, σ)
# shocks = rand(dist, h)
# shocks[:, 2:end] .= 0.
shocks = zeros(n_shocks_exogenous(m), h)
shocks[m.exogenous_shocks[:b_sh], 1] = -.05


_, ngdp_obs, _, _ = forecast(ngdp_sys, init_state, shocks)
_, cal_obs, _, _ = forecast(cal_sys, init_state, shocks)

using Plots
p3 = plot(1:h, ngdp_obs[m.observables[:obs_gdp], :], label = "Permanent NGDP", linewidth = 3)
plot!(1:h, cal_obs[m.observables[:obs_gdp], :], label = "Regime-Switching NGDP, P[NGDP] = 0.95",
      legend = :bottomright, linewidth = 3)

# Test 4
m <= Setting(:alternative_policy, DSGE.ngdp())
Tcal, Rcal, Ccal = gensys_prob_regime_switch(m, [0.05, 0.95], altpolicies)
Tcal, Rcal, Ccal = DSGE.augment_states(m, Tcal, Rcal, Ccal)
Tngdp, Rngdp, Cngdp = altpolicies[1].solve(m)
ngdp_sys = deepcopy(sys)
cal_sys  = deepcopy(sys)
ngdp_sys.transition.TTT[:, :] = Tngdp
ngdp_sys.transition.RRR[:, :] = Rngdp
ngdp_sys.transition.CCC[:, :] = Cngdp
cal_sys.transition.TTT[:, :] = Tcal
cal_sys.transition.RRR[:, :] = Rcal
cal_sys.transition.CCC[:, :] = Ccal
init_state = zeros(size(Tcal, 1))
h = 40
# Random.seed!(2020)
# μ = zeros(n_shocks_exogenous(m))
# σ = sqrt.(cal_sys[:QQ])
# dist = DegenerateMvNormal(μ, σ)
# shocks = rand(dist, h)
# shocks[:, 2:end] .= 0.
shocks = zeros(n_shocks_exogenous(m), h)
shocks[m.exogenous_shocks[:b_sh], 1] = -.05


_, ngdp_obs, _, _ = forecast(ngdp_sys, init_state, shocks)
_, cal_obs, _, _ = forecast(cal_sys, init_state, shocks)

using Plots
p4 = plot(1:h, ngdp_obs[m.observables[:obs_gdp], :], label = "Permanent NGDP", linewidth = 3)
plot!(1:h, cal_obs[m.observables[:obs_gdp], :], label = "Regime-Switching NGDP, P[NGDP] = 0.95",
      legend = :bottomright, linewidth = 3)
=#
