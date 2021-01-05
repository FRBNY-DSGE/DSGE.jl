"""
```
solve(m::AbstractDSGEModel; apply_altpolicy = false)
```

Driver to compute the model solution and augment transition matrices.

### Inputs

- `m`: the model object

## Keyword Arguments

- `apply_altpolicy::Bool`: whether or not to solve the model under the
  alternative policy. This should be `true` when we solve the model to
  forecast, but `false` when computing smoothed historical states (since
  the past was estimated under the baseline rule).
- `regime_switching::Bool`: true if the state space system features regime switching
- `regimes::Union{Int, Vector{Int}, UnitRange{Int}}`: specifies the specific regime to solve for.

### Outputs
 - TTT, RRR, and CCC matrices of the state transition equation:
```
S_t = TTT*S_{t-1} + RRR*ϵ_t + CCC
```
"""
function solve(m::AbstractDSGEModel{T}; apply_altpolicy = false,
               regime_switching::Bool = false,
               gensys_regimes::Vector{UnitRange{Int64}} = UnitRange{Int64}[1:1],
               gensys2_regimes::Vector{UnitRange{Int64}} = UnitRange{Int64}[1:1],
               regimes::Vector{Int} = Int[1],
               verbose::Symbol = :high) where {T <: Real}

    altpolicy_solve = alternative_policy(m).solve
    uncertain_altpolicy = haskey(get_settings(m), :uncertain_altpolicy) ? get_setting(m, :uncertain_altpolicy) : false

    if regime_switching
        return solve_regime_switching(m; apply_altpolicy = apply_altpolicy, gensys_regimes = gensys_regimes,
                                      gensys2_regimes = gensys2_regimes, regimes = regimes,
                                      uncertain_altpolicy = uncertain_altpolicy, verbose = verbose)
    else
        if get_setting(m, :solution_method) == :gensys
            if altpolicy_solve == solve || !apply_altpolicy

                # Get equilibrium condition matrices
                Γ0, Γ1, C, Ψ, Π  = eqcond(m)

                # Solve model
                TTT_gensys, CCC_gensys, RRR_gensys, eu = gensys(Γ0, Γ1, C, Ψ, Π, 1+1e-6, verbose = verbose)

                # Check for LAPACK exception, existence and uniqueness
                if eu[1] != 1 || eu[2] != 1
                    throw(GensysError())
                end

                TTT_gensys = real(TTT_gensys)
                RRR_gensys = real(RRR_gensys)
                CCC_gensys = real(CCC_gensys)

                # Augment states
                TTT, RRR, CCC = augment_states(m, TTT_gensys, RRR_gensys, CCC_gensys)
            else
                # Change the policy rule
                TTT, RRR, CCC = altpolicy_solve(m)
            end

            if uncertain_altpolicy && apply_altpolicy
                weights = get_setting(m, :alternative_policy_weights)
                altpols = get_setting(m, :alternative_policies)
                inds = 1:n_states(m)
                TTT_gensys, RRR_gensys, CCC_gensys = gensys_uncertain_altpol(m, weights, altpols; apply_altpolicy = apply_altpolicy,
                                                                             TTT = TTT[inds, inds])

                # Augment states
                TTT, RRR, CCC = augment_states(m, TTT_gensys, RRR_gensys, CCC_gensys)
            end
        elseif get_setting(m, :solution_method) == :klein
            TTT_jump, TTT_state = klein(m)

            # Transition
            TTT, RRR = klein_transition_matrices(m, TTT_state, TTT_jump)
            CCC = zeros(n_model_states(m))
        end

        return TTT, RRR, CCC
    end
end

"""
solve(m::PoolModel)
```

Driver to compute the model solution when using the PoolModel type

### Inputs

- `m`: the PoolModel object

### Outputs
- Φ: transition function
- F_ϵ: distribution of structural shock
- F_λ: prior on the initial λ_0

"""
function solve(m::PoolModel)
    return transition(m)
end

"""
```
GensysError <: Exception
```
A `GensysError` is thrown when:

1. Gensys does not give existence and uniqueness, or
2. A LAPACK error was thrown while computing the Schur decomposition of Γ0 and Γ1

If a `GensysError`is thrown during Metropolis-Hastings, it is caught by
`posterior`.  `posterior` then returns a value of `-Inf`, which
Metropolis-Hastings always rejects.

### Fields

* `msg::String`: Info message. Default = \"Error in Gensys\"
"""
mutable struct GensysError <: Exception
    msg::String
end
GensysError() = GensysError("Error in Gensys")
Base.showerror(io::IO, ex::GensysError) = print(io, ex.msg)

"""
```
KleinError <: Exception
```
A `KleinError` is thrown when:
1. A LAPACK error was thrown while computing the pseudo-inverse of U22*U22'

If a `KleinError`is thrown during Metropolis-Hastings, it is caught by
`posterior`.  `posterior` then returns a value of `-Inf`, which
Metropolis-Hastings always rejects.

### Fields

* `msg::String`: Info message. Default = \"Error in Klein\"
"""
mutable struct KleinError <: Exception
    msg::String
end
KleinError() = KleinError("Error in Klein")
Base.showerror(io::IO, ex::KleinError) = print(io, ex.msg)

"""
```
solve_regime_switching(m::AbstractDSGEModel; apply_altpolicy = false)

solve_one_regime(m::AbstractDSGEModel; apply_altpolicy = false)

solve_non_gensys2_regimes!(m::AbstractDSGEModel, Γ0s::Vector{Matrix{S}}, Γ01::Vector{Matrix{S}},
    Cs::Vector{Vector{S}}, Ψs::Vector{Matrix{S}}, Πs::Vector{Matrix{S}};
    TTTs::Vector{Matrix{S}}, RRRs::Vector{Matrix{S}}, CCCs::Vector{Vector{S}};
    regimes::Vector{Int} = Int[1], uncertain_altpolicy::Bool = false,
    apply_altpolicy::Bool = false, altpolicy_solve::Function = solve,
    verbose::Symbol = :high)

solve_gensys2!(m::AbstractDSGEModel, Γ0s::Vector{Matrix{S}}, Γ01::Vector{Matrix{S}},
    Cs::Vector{Vector{S}}, Ψs::Vector{Matrix{S}}, Πs::Vector{Matrix{S}},
    TTT_gensys_final::AbstractMatrix{S}, RRR_gensys_final::AbstractMatrix{S},
    CCC_gensys_final::AbstractMatrix{S},
    TTTs::Vector{Matrix{S}}, RRRs::Vector{Matrix{S}}, CCCs::Vector{Vector{S}};
    fcast_gensys2_regimes::Vector{Int} = Int[1], verbose::Symbol = :high)
```

calculates the reduced form state transition matrices in the case of regime switching.
These functions are intended to be internal functions hidden from the user but are separate
from the definition of the main `solve` function to ensure the main function is
comprehensible.
"""
function solve_regime_switching(m::AbstractDSGEModel{T}; apply_altpolicy::Bool = false,
                                uncertain_altpolicy::Bool = false,
                                gensys_regimes::Vector{UnitRange{Int64}} = UnitRange{Int64}[1:1],
                                gensys2_regimes::Vector{UnitRange{Int64}} = UnitRange{Int64}[1:1],
                                regimes::Vector{Int} = Int[1],
                                verbose::Symbol = :high) where {T <: Real}

    altpolicy_solve = alternative_policy(m).solve
    gensys2 = haskey(get_settings(m), :gensys2) ? get_setting(m, :gensys2) : false
    uncertain_zlb = haskey(get_settings(m), :uncertain_zlb) ? get_setting(m, :uncertain_zlb) : false
    if get_setting(m, :solution_method) == :gensys
        if length(regimes) == 1 # Calculate the solution to a specific regime
            return solve_one_regime(m; apply_altpolicy = apply_altpolicy, regime = regimes[1],
                                    uncertain_altpolicy = uncertain_altpolicy, verbose = verbose)
        else # Calculate the reduced-form state space matrices for all regimes
            Γ0s = Vector{Matrix{Float64}}(undef, length(regimes))
            Γ1s = Vector{Matrix{Float64}}(undef, length(regimes))
            Cs = Vector{Vector{Float64}}(undef, length(regimes))
            Ψs = Vector{Matrix{Float64}}(undef, length(regimes))
            Πs = Vector{Matrix{Float64}}(undef, length(regimes))

            for reg in regimes
                Γ0s[reg], Γ1s[reg], Cs[reg], Ψs[reg], Πs[reg] = eqcond(m, reg)
            end

            TTTs = Vector{Matrix{Float64}}(undef, length(regimes))
            RRRs = Vector{Matrix{Float64}}(undef, length(regimes))
            CCCs = Vector{Vector{Float64}}(undef, length(regimes))

            # Solve model for gensys regimes
            for reg_range in gensys_regimes
                solve_non_gensys2_regimes!(m, Γ0s, Γ1s, Cs, Ψs, Πs, TTTs, RRRs, CCCs;
                                           regimes = collect(reg_range),
                                           apply_altpolicy = false,
                                           verbose = verbose)
            end

            if gensys2
                for reg_range in gensys2_regimes
                    solve_gensys2!(m, Γ0s, Γ1s, Cs, Ψs, Πs,
                                   TTTs, RRRs, CCCs; gensys2_regimes = collect(reg_range), uncertain_zlb = uncertain_zlb,
                                   verbose = verbose)
                    # TODO: extend "uncertain ZLB" to "uncertain_gensys2" since it can in principle be used for
                    #       any temporary policy
                end
            end

            return TTTs, RRRs, CCCs
        end
    else
        error("Regime switching has not been implemented for other solution methods.")
    end
end

function solve_one_regime(m::AbstractDSGEModel{T}; apply_altpolicy = false,
                          regime::Int = 1, uncertain_altpolicy::Bool = false,
                          verbose::Symbol = :high) where {T <: Real}

    altpolicy_solve = alternative_policy(m).solve

    if altpolicy_solve == solve || !apply_altpolicy

        # Get equilibrium condition matrices
        Γ0, Γ1, C, Ψ, Π  = eqcond(m, regime)

        # Solve model
        TTT_gensys, CCC_gensys, RRR_gensys, eu = gensys(Γ0, Γ1, C, Ψ, Π, 1+1e-6, verbose = verbose)

        # Check for LAPACK exception, existence and uniqueness
        if eu[1] != 1 || eu[2] != 1
            throw(GensysError())
        end

        # Since we only have one regime in solve_one_regime, no time-varying credibility
        if (haskey(get_settings(m), :flexible_ait_policy_change) ?
            get_setting(m, :flexible_ait_policy_change) : false) &&
            get_setting(m, :regime_dates)[regime] >= get_setting(m, :flexible_ait_policy_change_date)
            TTT_gensys, RRR_gensys, CCC_gensys =
            gensys_uncertain_altpol(m, get_setting(m, :imperfect_awareness_weights),
                                    AltPolicy[get_setting(m, :imperfect_awareness_historical_policy)];
                                    apply_altpolicy = apply_altpolicy, TTT = TTT_gensys,
                                    regime_switching = true, regimes = Int[regime])
        end

        TTT_gensys = real(TTT_gensys)
        RRR_gensys = real(RRR_gensys)
        CCC_gensys = real(CCC_gensys)

        # Augment states
        TTT, RRR, CCC = augment_states(m, TTT_gensys, RRR_gensys, CCC_gensys; regime_switching = true,
                                       reg = regime)
    else
        TTT, RRR, CCC = altpolicy_solve(m; regime_switching = true, regimes = Int[regime])
    end

    if uncertain_altpolicy && apply_altpolicy
        # Since we only have one regime in solve_one_regime, no time-varying credibility
        weights = get_setting(m, :alternative_policy_weights)
        altpols = get_setting(m, :alternative_policies)

        inds = 1:n_states(m)
        TTT_gensys, RRR_gensys, CCC_gensys = gensys_uncertain_altpol(m, weights, altpols; apply_altpolicy = apply_altpolicy,
                                                                     TTT = TTT[inds, inds])

        # Augment states
        TTT, RRR, CCC = augment_states(m, TTT_gensys, RRR_gensys, CCC_gensys; regime_switching = true,
                                       regime = regime)
    end

    return TTT, RRR, CCC
end

function solve_non_gensys2_regimes!(m::AbstractDSGEModel, Γ0s::Vector{Matrix{S}}, Γ1s::Vector{Matrix{S}},
                                    Cs::Vector{Vector{S}}, Ψs::Vector{Matrix{S}}, Πs::Vector{Matrix{S}},
                                    TTTs::Vector{Matrix{S}}, RRRs::Vector{Matrix{S}}, CCCs::Vector{Vector{S}};
                                    regimes::Vector{Int} = Int[1],
                                    uncertain_altpolicy::Bool = false, apply_altpolicy::Bool = false,
                                    altpolicy_solve::Function = solve,
                                    verbose::Symbol = :high) where {S <: Real}

    for reg in regimes
        TTT_gensys, CCC_gensys, RRR_gensys, eu =
        gensys(Γ0s[reg], Γ1s[reg], Cs[reg], Ψs[reg], Πs[reg],
               1+1e-6, verbose = verbose)

        # Check for LAPACK exception, existence and uniqueness
        if eu[1] != 1 || eu[2] != 1
            throw(GensysError("Error in Gensys, Regime $reg"))
        end

        if get_setting(m, :gensys2) && get_setting(m, :uncertain_altpolicy)
            first_gensys2_regime = min(collect(keys(get_setting(m, :replace_eqcond_func_dict)))...)
            weights = get_setting(m, :imperfect_awareness_weights)
            histpol = AltPolicy[get_setting(m, :imperfect_awareness_historical_policy)]
            if reg >= first_gensys2_regime
                # Time-varying credibility weights for the regime reg

                if haskey(get_settings(m), :imperfect_awareness_varying_weights)
                    #weights = [get_setting(m, :imperfect_awareness_varying_weights)[i][reg-get_setting(m, :reg_forecast_start)+1]
                #               for i in 1:length(get_setting(m, :imperfect_awareness_varying_weights))]
                    #append!(weights, 1.0 - sum(weights))
                    weights = get_setting(m, :imperfect_awareness_varying_weights)[reg]
                end
                TTT_gensys, RRR_gensys, CCC_gensys =
                gensys_uncertain_altpol(m, weights, histpol; apply_altpolicy = true,
                                        TTT = TTT_gensys, regime_switching = true, regimes = Int[reg])
            end
        end

        TTT_gensys = real(TTT_gensys)
        RRR_gensys = real(RRR_gensys)
        CCC_gensys = real(CCC_gensys)

        # Populate the TTTs, etc., for regime `reg`
        TTTs[reg], RRRs[reg], CCCs[reg] =
        augment_states(m, TTT_gensys, RRR_gensys, CCC_gensys)
    end
    return TTTs, RRRs, CCCs
end

function solve_gensys2!(m::AbstractDSGEModel, Γ0s::Vector{Matrix{S}}, Γ1s::Vector{Matrix{S}},
                        Cs::Vector{Vector{S}}, Ψs::Vector{Matrix{S}}, Πs::Vector{Matrix{S}},
                        TTTs::Vector{Matrix{S}}, RRRs::Vector{Matrix{S}}, CCCs::Vector{Vector{S}};
                        gensys2_regimes::Vector{Int} = Int[1], uncertain_zlb::Bool = false,
                        verbose::Symbol = :high) where {S <: Real}
    # Solve for the final regime of the alternative rule
    altpolicy_solve = get_setting(m, :alternative_policy).solve
    TTT_final, RRR_final, CCC_final = altpolicy_solve(m; regime_switching = true,
                                                      regimes = Int[last(gensys2_regimes)])

    n_endo = length(m.endogenous_states)
    TTT_final = TTT_final[1:n_endo, 1:n_endo] # make sure the non-augmented version
    RRR_final = RRR_final[1:n_endo, :]        # is returned
    CCC_final = CCC_final[1:n_endo]

    # if we're using an uncertain alternative policy, we have to compute the
    # weighted final transition matrix
    uncertain_altpolicy = haskey(get_settings(m), :uncertain_altpolicy) ? get_setting(m, :uncertain_altpolicy) : false
    if uncertain_altpolicy
        if haskey(get_settings(m), :imperfect_awareness_varying_weights)
            weights = get_setting(m, :imperfect_awareness_varying_weights)
            #append!(weights, 1.0 - sum(weights))
            weights = weights[last(gensys2_regimes)]
        else
            weights = get_setting(m, :imperfect_awareness_weights)
        end

        altpols = get_setting(m, :alternative_policies)

        TTT_final_weighted, RRR_final_weighted, CCC_final_weighted =
        gensys_uncertain_altpol(m, weights, altpols; apply_altpolicy = true,
                                TTT = TTT_final, regimes = Int[last(gensys2_regimes)])

        TTT_final_weighted = real(TTT_final_weighted)
        RRR_final_weighted = real(RRR_final_weighted)
        CCC_final_weighted = real(CCC_final_weighted)
    end

    TTT_final = real(TTT_final)
    RRR_final = real(RRR_final)
    CCC_final = real(CCC_final)

    # Populate TTTs, RRRs, CCCs matrices
    if uncertain_zlb
        # Setup
        ffreg = first(gensys2_regimes) + 1
        altpols, weights = if haskey(get_settings(m), :imperfect_awareness_varying_weights) &&
            haskey(get_settings(m), :alternative_policies)
            get_setting(m, :alternative_policies), get_setting(m, :imperfect_awareness_varying_weights)
        elseif haskey(get_settings(m), :alternative_policies) && haskey(get_settings(m), :imperfect_awareness_weights)
            get_setting(m, :alternative_policies), get_setting(m, :imperfect_awareness_weights)
        else
            error("Alternative policies and/or weights were not specified.")
        end
        if isa(weights, AbstractDict)
            # populate_reg should work for now b/c should coincide w/entries of replace_eqcond_func_dict
            weights = [weights[i] for i in first(gensys2_regimes)+1:last(gensys2_regimes)]
        end

        @assert length(altpols) == 1 "Currently, uncertain_zlb works only for two policies (two possible MP rules)."
        Talt, _, Calt = altpols[1].solve(m)

        # Calculate gensys2 matrices under belief that the desired lift-off policy will occur
        # FIGURE OUT HOW TO NOT USE all gensys2_regimes, just use the minimal number of required eqcond matrices
        Tcal, Rcal, Ccal = gensys_cplus(m, Γ0s[gensys2_regimes], Γ1s[gensys2_regimes],
                                        Cs[gensys2_regimes], Ψs[gensys2_regimes], Πs[gensys2_regimes],
                                        TTT_final, RRR_final, CCC_final,
                                        T_switch = length(gensys2_regimes)-1)
        Tcal[end] = TTT_final
        Rcal[end] = RRR_final
        Ccal[end] = CCC_final

        # Now calculate transition matrices under an uncertain ZLB
        # Γ0_til, etc., are eqcond matrices implementing ZLB, i.e. zero_rate_rule
        # It is assumed that there is no re
        Γ0_til, Γ1_til, Γ2_til, C_til, Ψ_til =
        gensys_to_predictable_form(Γ0s[ffreg], Γ1s[ffreg], Cs[ffreg], Ψs[ffreg], Πs[ffreg])

        ng2  = length(Tcal) - 1 # number of gensys2 regimes
        nzlb = haskey(get_settings(m), :temporary_zlb_length) ? get_setting(m, :temporary_zlb_length) : ng2

        # Use Tcal, Rcal, & Ccal from 2 as inputs b/c use t + 1 matrix, not t
        # Then, if nzlb = 1, Tcal should have length 2, and you only need the lift-off matrix
        Tcal[1:(1 + nzlb)], Rcal[1:(1 + nzlb)], Ccal[1:(1 + nzlb)] =
        gensys_uncertain_zlb(weights, Talt[1:n_endo, 1:n_endo], Calt[1:n_endo],
                             Tcal[2:(1 + nzlb)], Rcal[2:(1 + nzlb)], Ccal[2:(1 + nzlb)],
                             Γ0_til, Γ1_til, Γ2_til, C_til, Ψ_til)

        if nzlb != ng2
            @show "Are we ending up in this block?"
            # TODO: generalize this code block
            # It is currently assumed that ffreg is the first regime w/ZLB (and cannot be another temporary altpolicy),
            # so ffreg + nzlb is the first regime w/out ZLB. Starting from nzlb + 1 ensures we populate the lift-off
            # regime correctly, as the lift-off is no longer the final regime, and won't be covered by
            # setting Tcal[end] = TTT_gensys_final
            # First UnitRange is for actual regime number, second for the index of Tcal, Rcal, & Ccal
            for (fcast_reg, ical) in zip((ffreg + nzlb):gensys2_regimes[end], (nzlb + 1):(ng2 + 1))
                @show fcast_reg, ical
                TTT_gensys, CCC_gensys, RRR_gensys, eu =
                gensys(Γ0s[fcast_reg], Γ1s[fcast_reg], Cs[fcast_reg], Ψs[fcast_reg], Πs[fcast_reg],
                       1+1e-6, verbose = verbose)

                # Check for LAPACK exception, existence and uniqueness
                if eu[1] != 1 || eu[2] != 1
                    throw(GensysError("Error in Gensys, Regime $fcast_reg"))
                end
                TTT_gensys = real(TTT_gensys)
                RRR_gensys = real(RRR_gensys)
                CCC_gensys = real(CCC_gensys)

                if haskey(get_settings(m), :flexible_ait_policy_change) ?
                    get_setting(m, :flexible_ait_policy_change) : false &&
                    get_setting(m, :regime_dates)[fcast_reg] >= get_setting(m, :flexible_ait_policy_change_date)

                    if haskey(get_settings(m), :imperfect_awareness_varying_weights)
                        imperfect_wt = get_setting(m, :imperfect_awareness_varying_weights)
                        weights = [imperfect_wt[i][fcast_reg - ffreg + 1] for i in 1:length(imperfect_wt)]
                        append!(weights, 1.0 - sum(weights))
                    else
                        weights = get_setting(m, :imperfect_awareness_weights)
                    end
                    histpol = AltPolicy[get_setting(m, :imperfect_awareness_historical_policy)]
                    TTT_gensys, RRR_gensys, CCC_gensys =
                    gensys_uncertain_altpol(m, weights, histpol;
                                            TTT = TTT_gensys, regime_switching = true, regimes = Int[fcast_reg])
                elseif haskey(get_settings(m), :uncertain_altpolicy) ? get_setting(m, :uncertain_altpolicy) : false
                    # THE INDEXING HERE ASSUMES TEMORARY ZLB STARTS IN THE PERIOD AFTER THE CONDITIONAL FORECAST
                    # HORIZON ENDS AND THAT THE ALTERNATIVE POLICY VARYING WEIGHTS VECTOR SPECIFIES WEIGHTS
                    # STARTING IN THAT SAME PERIOD (i.e. alternative_policy_varying_weights[i][1] is the first ZLB period)
                    if haskey(get_settings(m), :imperfect_awareness_varying_weights)
                        weights = get_setting(m, :imperfect_awareness_varying_weights)[fcast_reg]
                    else
                        weights = get_setting(m, :imperfect_awareness_weights)
                    end

                    altpols = get_setting(m, :alternative_policies)

                    TTT_gensys, RRR_gensys, CCC_gensys =
                    gensys_uncertain_altpol(m, weights, altpols; apply_altpolicy = true, TTT = TTT_gensys,
                                            regime_switching = true, regimes = Int[fcast_reg])
                end

                Tcal[ical] = TTT_gensys
                Rcal[ical] = RRR_gensys
                Ccal[ical] = CCC_gensys
            end
        end

        if uncertain_altpolicy
            Tcal[end] = TTT_final_weighted
            Rcal[end] = RRR_final_weighted
            Ccal[end] = CCC_final_weighted
        else
            Tcal[end] = TTT_final
            Rcal[end] = RRR_final
            Ccal[end] = CCC_final
        end

    else
        Tcal, Rcal, Ccal = gensys_cplus(m, Γ0s[gensys2_regimes], Γ1s[gensys2_regimes],
                                        Cs[gensys2_regimes], Ψs[gensys2_regimes], Πs[gensys2_regimes],
                                        TTT_final, RRR_final, CCC_final;
                                        T_switch = length(gensys2_regimes)-1)

        if uncertain_altpolicy
            Tcal[end] = TTT_final_weighted
            Rcal[end] = RRR_final_weighted
            Ccal[end] = CCC_final_weighted
        else
            Tcal[end] = TTT_final
            Rcal[end] = RRR_final
            Ccal[end] = CCC_final
        end
    end

    # only need to populate regimes during which a temporary altpolicy holds
    populate_reg = gensys2_regimes[2]:gensys2_regimes[end]
    @show populate_reg
    for (i, reg) in enumerate(populate_reg)
        TTTs[reg], RRRs[reg], CCCs[reg] = augment_states(m, Tcal[i], Rcal[i], Ccal[i])
    end

    return TTTs, RRRs, CCCs
end
