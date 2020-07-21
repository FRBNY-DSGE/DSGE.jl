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
               hist_regimes::Vector{Int} = Int[1],
               fcast_regimes::Vector{Int} = Int[1],
               regimes::Vector{Int} = Int[1],
               verbose::Symbol = :high) where {T <: Real}

    altpolicy_solve = alternative_policy(m).solve
    uncertain_altpolicy = haskey(get_settings(m), :uncertain_altpolicy) ? get_setting(m, :uncertain_altpolicy) : false

    if regime_switching
        return solve_regime_switching(m; apply_altpolicy = apply_altpolicy, hist_regimes = hist_regimes,
                                      fcast_regimes = fcast_regimes, regimes = regimes, uncertain_altpolicy = uncertain_altpolicy,
                                      verbose = verbose)
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

solve_histregimes!(m::AbstractDSGEModel, Γ0s::Vector{Matrix{S}}, Γ01::Vector{Matrix{S}},
    Cs::Vector{Vector{S}}, Ψs::Vector{Matrix{S}}, Πs::Vector{Matrix{S}};
    TTTs::Vector{Matrix{S}}, RRRs::Vector{Matrix{S}}, CCCs::Vector{Vector{S}};
    n_hist_regimes::Vector{Int} = Int[1], verbose::Symbol = :high)

solve_fcastregimes!(m::AbstractDSGEModel, Γ0s::Vector{Matrix{S}}, Γ01::Vector{Matrix{S}},
    Cs::Vector{Vector{S}}, Ψs::Vector{Matrix{S}}, Πs::Vector{Matrix{S}};
    TTTs::Vector{Matrix{S}}, RRRs::Vector{Matrix{S}}, CCCs::Vector{Vector{S}};
    n_fcast_regimes::Vector{Int} = Int[1], apply_altpolicy::Bool = false,
    altpolicy_solve = altpolicy_solve, verbose::Symbol = :high)

solve_gensys2!(m::AbstractDSGEModel, Γ0s::Vector{Matrix{S}}, Γ01::Vector{Matrix{S}},
    Cs::Vector{Vector{S}}, Ψs::Vector{Matrix{S}}, Πs::Vector{Matrix{S}},
    TTT_gensys_final::AbstractMatrix{S}, RRR_gensys_final::AbstractMatrix{S},
    CCC_gensys_final::AbstractMatrix{S},
    TTTs::Vector{Matrix{S}}, RRRs::Vector{Matrix{S}}, CCCs::Vector{Vector{S}};
    fcast_regimes::Vector{Int} = Int[1], verbose::Symbol = :high)
```

calculates the reduced form state transition matrices in the case of regime switching.
These functions are intended to be internal functions hidden from the user but are separate
from the definition of the main `solve` function to ensure the main function is
comprehensible.
"""
function solve_regime_switching(m::AbstractDSGEModel{T}; apply_altpolicy::Bool = false,
                                uncertain_altpolicy::Bool = false,
                                hist_regimes::Vector{Int} = Int[1],
                                fcast_regimes::Vector{Int} = Int[1],
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

            # Solve model for regimes corresponding to the history of a forecast
            solve_histregimes!(m, Γ0s, Γ1s, Cs, Ψs, Πs, TTTs, RRRs, CCCs;
                               hist_regimes = hist_regimes, verbose = verbose)

            # Get the last state space matrices one period after the alternative policy ends (if there is one)
            if altpolicy_solve == solve || !apply_altpolicy
                # If normal rule
                TTT_gensys_final, CCC_gensys_final, RRR_gensys_final, eu = gensys(Γ0s[end], Γ1s[end], Cs[end],
                                                                                  Ψs[end], Πs[end],
                                                                                  1+1e-6, verbose = verbose)
                # Check for LAPACK exception, existence and uniqueness
                if eu[1] != 1 || eu[2] != 1
                    throw(GensysError("Error in Gensys, Regime $reg"))
                end
            else
                # If alternative rule
                n_endo = length(keys(m.endogenous_states))
                TTT_gensys_final, RRR_gensys_final, CCC_gensys_final = altpolicy_solve(m; regime_switching = true,
                                                                                       regimes = Int[get_setting(m, :n_regimes)])
                TTT_gensys_final = TTT_gensys_final[1:n_endo, 1:n_endo] # make sure the non-augmented version
                RRR_gensys_final = RRR_gensys_final[1:n_endo, :]        # is returned
                CCC_gensys_final = CCC_gensys_final[1:n_endo]
            end

            if uncertain_altpolicy && apply_altpolicy
                weights = get_setting(m, :alternative_policy_weights)
                altpols = get_setting(m, :alternative_policies)

                TTT_gensys_final, RRR_gensys_final, CCC_gensys_final =
                    gensys_uncertain_altpol(m, weights, altpols; apply_altpolicy = apply_altpolicy,
                                            TTT = TTT_gensys_final)
            end


            TTT_gensys_final = real(TTT_gensys_final)
            RRR_gensys_final = real(RRR_gensys_final)
            CCC_gensys_final = real(CCC_gensys_final)

            # Are there temporary policies in the forecast that had been unanticipated in the history?
            if gensys2
                # TODO: generalize this policy beyond assuming the final regime HAS to be the final regime-switch
                # and that the final regime is the policy after the temporary alternative policy finishes
                solve_gensys2!(m, Γ0s, Γ1s, Cs, Ψs, Πs, TTT_gensys_final, RRR_gensys_final, CCC_gensys_final,
                               TTTs, RRRs, CCCs; fcast_regimes = fcast_regimes, uncertain_zlb = uncertain_zlb,
                               verbose = verbose)
                # TODO: extend "uncertain ZLB" to "uncertain_gensys2" since it can in principle be used for
                #       any temporary policy
            else
                solve_fcastregimes!(m, Γ0s, Γ1s, Cs, Ψs, Πs, TTTs, RRRs, CCCs;
                                    fcast_regimes = fcast_regimes, apply_altpolicy = apply_altpolicy,
                                    altpolicy_solve = altpolicy_solve, verbose = verbose)
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

function solve_histregimes!(m::AbstractDSGEModel, Γ0s::Vector{Matrix{S}}, Γ1s::Vector{Matrix{S}},
                            Cs::Vector{Vector{S}}, Ψs::Vector{Matrix{S}}, Πs::Vector{Matrix{S}},
                            TTTs::Vector{Matrix{S}}, RRRs::Vector{Matrix{S}}, CCCs::Vector{Vector{S}};
                            hist_regimes::Vector{Int} = Int[1], verbose::Symbol = :high) where {S <: Real}
    for hist_reg in hist_regimes
        TTT_gensys, CCC_gensys, RRR_gensys, eu =
            gensys(Γ0s[hist_reg], Γ1s[hist_reg], Cs[hist_reg], Ψs[hist_reg], Πs[hist_reg],
                   1+1e-6, verbose = verbose)

        # Check for LAPACK exception, existence and uniqueness
        if eu[1] != 1 || eu[2] != 1
            throw(GensysError("Error in Gensys, Regime $hist_reg"))
        end
        TTT_gensys = real(TTT_gensys)
        RRR_gensys = real(RRR_gensys)
        CCC_gensys = real(CCC_gensys)

        # Put the usual system in for the history (smoothing)
        TTTs[hist_reg], RRRs[hist_reg], CCCs[hist_reg] =
            augment_states(m, TTT_gensys, RRR_gensys, CCC_gensys)
    end

    return TTTs, RRRs, CCCs
end

function solve_fcastregimes!(m::AbstractDSGEModel, Γ0s::Vector{Matrix{S}}, Γ1s::Vector{Matrix{S}},
                             Cs::Vector{Vector{S}}, Ψs::Vector{Matrix{S}}, Πs::Vector{Matrix{S}},
                             TTTs::Vector{Matrix{S}}, RRRs::Vector{Matrix{S}}, CCCs::Vector{Vector{S}};
                             fcast_regimes::Vector{Int} = Int[1], apply_altpolicy::Bool = false,
                             verbose::Symbol = :high, altpolicy_solve = nothing) where {S <: Real}

    for fcast_reg in fcast_regimes
        if altpolicy_solve == solve || !apply_altpolicy
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

            TTTs[fcast_reg], RRRs[fcast_reg], CCCs[fcast_reg] =
                augment_states(m, TTT_gensys, RRR_gensys, CCC_gensys)
        else
            TTTs[fcast_reg], RRRs[fcast_reg], CCCs[fcast_reg] =
                altpolicy_solve(m; regime_switching = true, regimes = Int[fcast_reg])
        end
    end

    return TTTs, RRRs, CCCs
end

function solve_gensys2!(m::AbstractDSGEModel, Γ0s::Vector{Matrix{S}}, Γ1s::Vector{Matrix{S}},
                        Cs::Vector{Vector{S}}, Ψs::Vector{Matrix{S}}, Πs::Vector{Matrix{S}},
                        TTT_gensys_final::AbstractMatrix{S}, RRR_gensys_final::AbstractMatrix{S},
                        CCC_gensys_final::AbstractVector{S},
                        TTTs::Vector{Matrix{S}}, RRRs::Vector{Matrix{S}}, CCCs::Vector{Vector{S}};
                        fcast_regimes::Vector{Int} = Int[1], uncertain_zlb::Bool = false,
                        verbose::Symbol = :high) where {S <: Real}

    # Calculate the matrices for the temporary alternative policies
    gensys2_regimes = (first(fcast_regimes) - 1):last(fcast_regimes) # TODO: generalize to multiple times in which we need to impose temporary alternative policies

    # Are there periods between the first forecast period and first period with temporary alt policies?
    # For example, are there conditional periods?
    n_no_alt_reg = get_setting(m, :n_fcast_regimes) - get_setting(m, :n_rule_periods) - 1
    if n_no_alt_reg > 0
        # Get the T, R, C matrices between the first forecast period & first rule period
        for fcast_reg in first(fcast_regimes):(first(fcast_regimes) + n_no_alt_reg - 1)
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

            TTTs[fcast_reg], RRRs[fcast_reg], CCCs[fcast_reg] = DSGE.augment_states(m, TTT_gensys,
                                                                                    RRR_gensys, CCC_gensys)
        end
        populate_reg = (fcast_regimes[end] - get_setting(m, :n_rule_periods)):fcast_regimes[end]
    else
        populate_reg = fcast_regimes
    end

    # Populate TTTs, RRRs, CCCs matrices
    if uncertain_zlb
        # Setup
        ffreg = first(fcast_regimes) + n_no_alt_reg
        altpols = get_setting(m, :alternative_policies)
        weights = get_setting(m, :alternative_policy_weights)
        @assert length(altpols) == 1 "Currently, uncertain_zlb works only for two policies (two possible MP rules)."
        Talt, Ralt, Calt = altpols[1].solve(m)

        # Calculate the desired lift-off policy
        altpolicy_solve = get_setting(m, :alternative_policy).solve
        TTT_liftoff, RRR_liftoff, CCC_liftoff = altpolicy_solve(m; regime_switching = true,
                                                                regimes = Int[get_setting(m, :n_regimes)])

        n_endo = length(m.endogenous_states)
        TTT_liftoff = TTT_liftoff[1:n_endo, 1:n_endo] # make sure the non-augmented version
        RRR_liftoff = RRR_liftoff[1:n_endo, :]        # is returned
        CCC_liftoff = CCC_liftoff[1:n_endo]

        # Calculate gensys2 matrices under belief that the desired lift-off policy will occur
        Tcal, Rcal, Ccal = gensys_cplus(m, Γ0s[gensys2_regimes], Γ1s[gensys2_regimes],
                                        Cs[gensys2_regimes], Ψs[gensys2_regimes], Πs[gensys2_regimes],
                                        TTT_liftoff, RRR_liftoff, CCC_liftoff)
        Tcal[end] = TTT_liftoff
        Rcal[end] = RRR_liftoff
        Ccal[end] = CCC_liftoff

        # Now calculate transition matrices under an uncertain ZLB
        Γ0_til, Γ1_til, Γ2_til, C_til, Ψ_til =
            gensys_to_predictable_form(Γ0s[ffreg], Γ1s[ffreg], Cs[ffreg], Ψs[ffreg], Πs[ffreg])

        Tcal, Rcal, Ccal =
            gensys_uncertain_zlb(weights, Talt[1:n_endo, 1:n_endo], Calt[1:n_endo], Tcal[2:end], Rcal[2:end], Ccal[2:end],
                                 Γ0_til, Γ1_til, Γ2_til, C_til, Ψ_til)

        Tcal[end] = TTT_gensys_final
        Rcal[end] = RRR_gensys_final
        Ccal[end] = CCC_gensys_final

        for (i, reg) in enumerate(populate_reg)
            TTTs[reg], RRRs[reg], CCCs[reg] = augment_states(m, Tcal[i], Rcal[i], Ccal[i])
        end
    else
        Tcal, Rcal, Ccal = gensys_cplus(m, Γ0s[gensys2_regimes], Γ1s[gensys2_regimes],
                                        Cs[gensys2_regimes], Ψs[gensys2_regimes], Πs[gensys2_regimes],
                                        TTT_gensys_final, RRR_gensys_final, CCC_gensys_final)
        Tcal[end] = TTT_gensys_final
        Rcal[end] = RRR_gensys_final
        Ccal[end] = CCC_gensys_final

        for (i, fcast_reg) in enumerate(populate_reg)
            TTTs[fcast_reg], RRRs[fcast_reg], CCCs[fcast_reg] = augment_states(m, Tcal[i], Rcal[i], Ccal[i])
        end
    end

    return TTTs, RRRs, CCCs
end
