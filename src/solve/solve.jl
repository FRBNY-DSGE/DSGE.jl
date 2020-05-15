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
               regime_switching::Bool = false, regimes::Union{Int, Vector{Int}, UnitRange{Int}} = 1,
               verbose::Symbol = :high) where {T <: Real}
    altpolicy_solve = alternative_policy(m).solve

    # Need to skip that block of code (but this is clunky and confusing and need to fix)
    if get_setting(m, :gensys2)
        regime_switching = false
    end

    if regime_switching
        if get_setting(m, :solution_method) == :gensys
            if altpolicy_solve == solve || !apply_altpolicy
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
                    Γ0s = Vector{Matrix{T}}(undef, length(regimes))
                    Γ1s = Vector{Matrix{T}}(undef, length(regimes))
                    Cs = Vector{Vector{T}}(undef, length(regimes))
                    Ψs = Vector{Matrix{T}}(undef, length(regimes))
                    Πs = Vector{Matrix{T}}(undef, length(regimes))
                    for reg in regimes
                        Γ0s[reg], Γ1s[reg], Cs[reg], Ψs[reg], Πs[reg] = eqcond(m, reg)
                    end

                    n_regimes = length(regimes)
                    TTTs = Vector{Matrix{T}}(undef, n_regimes)
                    RRRs = Vector{Matrix{T}}(undef, n_regimes)
                    CCCs = Vector{Vector{T}}(undef, n_regimes)

                    # Solve model
                    for reg in regimes
                        TTT_gensys, CCC_gensys, RRR_gensys, eu =
                            gensys(Γ0s[reg], Γ1s[reg], Cs[reg], Ψs[reg], Πs[reg], 1+1e-6, verbose = verbose)

                        # Check for LAPACK exception, existence and uniqueness
                        if eu[1] != 1 || eu[2] != 1
                            throw(GensysError("Error in Gensys, Regime $reg"))
                        end

                        TTT_gensys = real(TTT_gensys)
                        RRR_gensys = real(RRR_gensys)
                        CCC_gensys = real(CCC_gensys)

                        # Augment states
                        TTTs[reg], RRRs[reg], CCCs[reg] = augment_states(m, TTT_gensys, RRR_gensys, CCC_gensys;
                                                                         regime_switching = regime_switching,
                        reg = reg)
                    end

                    return TTTs, RRRs, CCCs
                end
            else
                # Change the policy rule
                TTTs, RRRs, CCCs = altpolicy_solve(m; regime_switching = regime_switching,
                                                   regimes = regimes)
            end
        else
            # Change the policy rule
            TTT, RRR, CCC = altpolicy_solve(m)
        end
    elseif get_setting(m, :solution_method) == :klein
        TTT_jump, TTT_state = klein(m)

                if get_setting(m, :gensys2)
                    TTTs = Vector{Matrix{Float64}}(undef, get_setting(m, :n_regimes))
                    RRRs = Vector{Matrix{Float64}}(undef, get_setting(m, :n_regimes))
                    CCCs = Vector{Vector{Float64}}(undef, get_setting(m, :n_regimes))
                    # Put teh usual system in for the history (smoothing)
                    TTTs[1], RRRs[1], CCCs[1] = augment_states(m, TTT_gensys, RRR_gensys, CCC_gensys)

                    Γ0, Γ1, C, Ψ, Π  = eqcond(m, new_policy = true)
                    Tcal, Rcal, Ccal = gensys_cplus(m, Γ0, Γ1, C, Ψ, Π,
                                                    TTT_gensys, RRR_gensys, CCC_gensys)

                    for i in 1:length(Tcal)
                        TTTs[i+1], RRRs[i+1], CCCs[i+1] = augment_states(m, Tcal[i], Rcal[i], Ccal[i]) #=;
                                                                         regime_switching = regime_switching,
                                                                         reg = reg) =#
                    end
                    return TTTs, RRRs, CCCs
                end

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
