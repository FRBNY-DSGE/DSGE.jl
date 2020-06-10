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
               hist_regimes::Union{Int, Vector{Int}, UnitRange{Int}} = 1,
               fcast_regimes::Union{Int, Vector{Int}, UnitRange{Int}} = 1,
               regimes::Union{Int, Vector{Int}, UnitRange{Int}} = 1,
               verbose::Symbol = :high) where {T <: Real}
    altpolicy_solve = alternative_policy(m).solve

    gensys2 = haskey(get_settings(m), :gensys2) ? get_setting(m, :gensys2) : false

    if regime_switching
        if get_setting(m, :solution_method) == :gensys
            #if altpolicy_solve == solve || !apply_altpolicy
                if isa(regimes, Int)
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
                    TTT, RRR, CCC = augment_states(m, TTT_gensys, RRR_gensys, CCC_gensys; regime_switching = regime_switching,
                                                   regime = regimes)

                    return TTT, RRR, CCC
                else
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

                    # Solve model
                    for hist_reg in hist_regimes
                        TTT_gensys, CCC_gensys, RRR_gensys, eu =
                            gensys(Γ0s[hist_reg], Γ1s[hist_reg], Cs[hist_reg], Ψs[hist_reg], Πs[hist_reg],
                                   1+1e-6, verbose = verbose)

                        # Check for LAPACK exception, existence and uniqueness
                        if eu[1] != 1 || eu[2] != 1
                            throw(GensysError("Error in Gensys, Regime $reg"))
                        end

                        TTT_gensys = real(TTT_gensys)
                        RRR_gensys = real(RRR_gensys)
                        CCC_gensys = real(CCC_gensys)

                        # Put teh usual system in for the history (smoothing)
                        TTTs[hist_reg], RRRs[hist_reg], CCCs[hist_reg] = DSGE.augment_states(m, TTT_gensys,
                                                                                             RRR_gensys, CCC_gensys)
                    end

                    # Get the last state space matrices one period after the policy ends
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
                        TTT_gensys_final, RRR_gensys_final, CCC_gensys_final = altpolicy_solve(m; regime_switching = regime_switching,
                                                   regimes = get_setting(m, :n_regimes))
                        TTT_gensys_final = TTT_gensys_final[1:n_endo, 1:n_endo]
                        RRR_gensys_final = RRR_gensys_final[1:n_endo, :]
                        CCC_gensys_final = CCC_gensys_final[1:n_endo]
                    end


                    TTT_gensys_final = real(TTT_gensys_final)
                    RRR_gensys_final = real(RRR_gensys_final)
                    CCC_gensys_final = real(CCC_gensys_final)


                    for fcast_reg in fcast_regimes
                        Γ0s[fcast_reg], Γ1s[fcast_reg], Cs[fcast_reg], Ψs[fcast_reg], Πs[fcast_reg]  =
                            eqcond(m, fcast_reg, new_policy = gensys2)
                    end
                    if gensys2
                        gensys2_regimes = first(fcast_regimes)-1:last(fcast_regimes)
                        Tcal, Rcal, Ccal = DSGE.gensys_cplus(m, Γ0s[gensys2_regimes], Γ1s[gensys2_regimes],
                                                             Cs[gensys2_regimes], Ψs[gensys2_regimes], Πs[gensys2_regimes],
                                                             TTT_gensys_final,
                                                             RRR_gensys_final,
                                                             CCC_gensys_final)
                        Tcal[end] = TTT_gensys_final
                        Rcal[end] = RRR_gensys_final
                        Ccal[end] = CCC_gensys_final

                        for (i, fcast_reg) = enumerate(fcast_regimes)
                            TTTs[fcast_reg], RRRs[fcast_reg], CCCs[fcast_reg] = augment_states(m, Tcal[i], Rcal[i], Ccal[i])
                        end
                    else
                        for fcast_reg in fcast_regimes
                            if altpolicy_solve == solve || !apply_altpolicy
                                TTT_gensys, CCC_gensys, RRR_gensys, eu =
                                    gensys(Γ0s[fcast_reg], Γ1s[fcast_reg], Cs[fcast_reg], Ψs[fcast_reg], Πs[fcast_reg], 1+1e-6, verbose = verbose)

                                # Check for LAPACK exception, existence and uniqueness
                                if eu[1] != 1 || eu[2] != 1
                                    throw(GensysError("Error in Gensys, Regime $fcast_reg"))
                                end

                                TTT_gensys = real(TTT_gensys)
                                RRR_gensys = real(RRR_gensys)
                                CCC_gensys = real(CCC_gensys)

                                TTTs[fcast_reg], RRRs[fcast_reg], CCCs[fcast_reg] = augment_states(m, TTT_gensys, RRR_gensys, CCC_gensys) #=;
                                regime_switching = regime_switching,
                                reg = fcast_reg) =#
                            else
                                TTTs[fcast_reg], RRRs[fcast_reg], CCCs[fcast_reg] = altpolicy_solve(m; regime_switching = regime_switching,
                                                                  regimes = fcast_reg)
                            end
                        end
                    end
                    return TTTs, RRRs, CCCs
                end
          #=  else
                # Change the policy rule
                TTTs, RRRs, CCCs = altpolicy_solve(m; regime_switching = regime_switching,
                                                   regimes = regimes)
            end =#
        else
            # Change the policy rule
            TTT, RRR, CCC = altpolicy_solve(m)
        end
    elseif get_setting(m, :solution_method) == :klein
        TTT_jump, TTT_state = klein(m)

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

    return TTT, RRR, CCC
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
