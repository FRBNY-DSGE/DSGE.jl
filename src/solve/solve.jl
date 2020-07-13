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
               uncertain_altpolicy::Bool = false,
               verbose::Symbol = :high) where {T <: Real}

    altpolicy_solve = alternative_policy(m).solve

    if regime_switching
        return solve_regime_switching(m; apply_altpolicy = apply_altpolicy, hist_regimes = hist_regimes,
                                      fcast_regimes = fcast_regimes, regimes = regimes,
                                      uncertain_altpolicy = uncertain_altpolicy, verbose = verbose)
    else
        if get_setting(m, :solution_method) == :gensys
            if uncertain_altpolicy
                weights = get_setting(m, :alternative_policy_weights)
                altpols = get_setting(m, :alternative_policies)

                TTT_gensys, RRR_gensys, CCC_gensys = gensys_uncertain_altpol(m, weights, altpols; apply_altpolicy = apply_altpolicy)

                # Augment states
                TTT, RRR, CCC = augment_states(m, TTT_gensys, RRR_gensys, CCC_gensys)

            elseif altpolicy_solve == solve || !apply_altpolicy

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
        elseif get_setting(m, :solution_method) == :klein
            TTT_jump, TTT_state = klein(m)

            # Transition
            TTT, RRR = klein_transition_matrices(m, TTT_state, TTT_jump)
            CCC = zeros(n_model_states(m))
        end

        return TTT, RRR, CCC
    end
end

function uncertain_altpolicy_solve()

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
```

calculates the reduced form state transition matrices in the case of regime switching.
These functions are intended to be internal functions hidden from the user but are separate
from the definition of the main `solve` function to ensure the main function is
comprehensible.
"""
function solve_regime_switching(m::AbstractDSGEModel{T}; apply_altpolicy = false,
                                hist_regimes::Vector{Int} = Int[1],
                                fcast_regimes::Vector{Int} = Int[1],
                                regimes::Vector{Int} = Int[1],
                                uncertain_altpolicy::Bool = false,
                                verbose::Symbol = :high) where {T <: Real}

    altpolicy_solve = alternative_policy(m).solve
    gensys2 = haskey(get_settings(m), :gensys2) ? get_setting(m, :gensys2) : false

    if get_setting(m, :solution_method) == :gensys
        if length(regimes) == 1 # Calculate the solution to a specific regime
            solve_one_regime(m; apply_altpolicy = apply_altpolicy, regime = regimes[1],
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
                TTTs[hist_reg], RRRs[hist_reg], CCCs[hist_reg] = DSGE.augment_states(m, TTT_gensys,
                                                                                     RRR_gensys, CCC_gensys)
            end

            # Get the last state space matrices one period after the alternative policy ends (if there is one)
            if uncertain_altpolicy
                weights = get_setting(m, :alternative_policy_weights)
                altpols = get_setting(m, :alternative_policies)

                TTT_gensys_final, RRR_gensys_final, CCC_gensys_final =
                    gensys_uncertain_altpol(m, weights, altpols; apply_altpolicy = apply_altpolicy)

            elseif altpolicy_solve == solve || !apply_altpolicy
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
                TTT_gensys_final = TTT_gensys_final[1:n_endo, 1:n_endo]
                RRR_gensys_final = RRR_gensys_final[1:n_endo, :]
                CCC_gensys_final = CCC_gensys_final[1:n_endo]
            end

            TTT_gensys_final = real(TTT_gensys_final)
            RRR_gensys_final = real(RRR_gensys_final)
            CCC_gensys_final = real(CCC_gensys_final)

            # Now we handle regime switching during forecast time periods
            # 123-126 might be unnecessary now
            for fcast_reg in fcast_regimes
                Γ0s[fcast_reg], Γ1s[fcast_reg], Cs[fcast_reg], Ψs[fcast_reg], Πs[fcast_reg] =
                    eqcond(m, fcast_reg, new_policy = gensys2)
            end

            # Are there temporary policies in the forecast that had been unanticipated in the history?
            if gensys2
                # TODO: generalize this policy beyond assuming the final regime HAS to be the final regime-switch
                # and that the final regime is the policy after the temporary alternative policy finishes
                gensys2_regimes = (first(fcast_regimes) - 1):last(fcast_regimes) # TODO: generalize to multiple times in which we need to impose temporary alternative policies
                Tcal, Rcal, Ccal = gensys_cplus(m, Γ0s[gensys2_regimes], Γ1s[gensys2_regimes],
                                                Cs[gensys2_regimes], Ψs[gensys2_regimes], Πs[gensys2_regimes],
                                                TTT_gensys_final, RRR_gensys_final, CCC_gensys_final)
                Tcal[end] = TTT_gensys_final
                Rcal[end] = RRR_gensys_final
                Ccal[end] = CCC_gensys_final

                n_no_alt_reg = get_setting(m, :n_fcast_regimes) - get_setting(m, :n_rule_periods) - 1
                if n_no_alt_reg > 0
                    # Get the T, R, C matrices for the rule period + lift-off period
                    for (i, fcast_reg) in enumerate((gensys2_regimes[end] - # minus n_rule_periods b/c including lift-off regime
                                                     get_setting(m, :n_rule_periods)):gensys2_regimes[end])
                        TTTs[fcast_reg], RRRs[fcast_reg], CCCs[fcast_reg] = augment_states(m, Tcal[i], Rcal[i], Ccal[i])
                    end

                    # Get the T, R, C matrices between the first forecast period & first rule period
                    for fcast_reg in first(fcast_regimes):(first(fcast_regimes) + n_no_alt_reg - 1)
                        # TODO: place this code block inside its own function
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
                else
                    # The alternative policy will cover all the regimes in the forecast period
                    for (i, fcast_reg) in enumerate(fcast_regimes)
                        TTTs[fcast_reg], RRRs[fcast_reg], CCCs[fcast_reg] = augment_states(m, Tcal[i], Rcal[i], Ccal[i])
                    end
                end

                #=
                elseif uncertain_altpolicy
                # Regime-switching with the uncertain altpolicy approach
                TTTs[fcast_regimes], RRRs[fcast_regimes], CCCs[fcast_regimes] =
                gensys_uncertain_altpol(m, weights, altpols; apply_altpolicy = apply_altpolicy,
                regime_switching = regime_switching, regimes = fcast_regimes,
                Γ0s = Γ0s, Γ1s = Γ1s, Cs = Cs, Ψs = Ψs, Πs = Πs)
                for fcast_reg in fcast_regimes
                TTTs[fcast_reg], RRRs[fcast_reg], CCCs[fcast_reg] =
                augment_states(m, TTTs[fcast_reg], RRRs[fcast_reg], CCCs[fcast_reg])
                end
                =#
            else
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
            end
        end

        return TTTs, RRRs, CCCs
    else
        error("Regime switching has not been implemented for other solution methods.")
    end
end

function solve_one_regime(m::AbstractDSGEModel{T}; apply_altpolicy = false,
                          regime::Int = 1, uncertain_altpolicy::Bool = false,
                          verbose::Symbol = :high) where {T <: Real}

    altpolicy_solve = alternative_policy(m).solve

    if uncertain_altpolicy
        weights = get_setting(m, :alternative_policy_weights)
        altpols = get_setting(m, :alternative_policies)

        TTT_gensys, RRR_gensys, CCC_gensys = gensys_uncertain_altpol(m, weights, altpols; apply_altpolicy = apply_altpolicy)

        # Augment states
        TTT, RRR, CCC = augment_states(m, TTT_gensys, RRR_gensys, CCC_gensys; regime_switching = true,
                                       regime = regimes)
    elseif altpolicy_solve == solve || !apply_altpolicy

        # Get equilibrium condition matrices
        Γ0, Γ1, C, Ψ, Π  = eqcond(m, regimes)

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
                                       regime = regimes)

    else
        TTT, RRR, CCC = altpolicy_solve(m; regime_switching = true, regimes = regimes)
    end

    return TTT, RRR, CCC
end
