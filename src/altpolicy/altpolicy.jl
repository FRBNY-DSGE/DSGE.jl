"""
```
mutable struct AltPolicy
```

Type defining an alternative policy rule.

### Fields

- `key::Symbol`: alternative policy identifier

- `eqcond::Function`: a version of `DSGE.eqcond` which computes the equilibrium
  condition matrices under the alternative policy. Like `DSGE.eqcond`, it should
  take in one argument of mutable struct `AbstractDSGEModel` and return the `Γ0`,
  `Γ1`, `C`, `Ψ`, and `Π` matrices.

- `solve::Function`: a version of `DSGE.solve` which solves the model under the
  alternative policy. Like `DSGE.solve`, it should take in one argument of mutable
  struct `AbstractDSGEModel` and return the `TTT`, `RRR`, and `CCC` matrices.

- `replace_eqcond::Function`: Changes the equilibrium condition matrices to align
  with the given alternative policy.

- `forecast_init::Function`: a function that initializes forecasts under the
  alternative policy rule. Specifically, it accepts a model, an `nshocks` x
  `n_forecast_periods` matrix of shocks to be applied in the forecast, and a
  vector of initial states for the forecast. It must return a new matrix of
  shocks and a new initial state vector. If no adjustments to shocks or initial
  state vectors are necessary under the policy rule, this field may be omitted.

- `color::Colorant`: color to plot this alternative policy in. Defaults to blue.

- `linestyle::Symbol`: line style for forecast plots under this alternative
  policy. See options from `Plots.jl`. Defaults to `:solid`.

"""
mutable struct AltPolicy
    key::Symbol
    eqcond::Function
    solve::Function
    replace_eqcond::Function
    forecast_init::Function
    color::Colorant
    linestyle::Symbol
end

function AltPolicy(key::Symbol, eqcond_fcn::Function, solve::Function;
                   replace_eqcond::Function = identity,
                   forecast_init::Function = identity,
                   color::Colorant = RGB(0., 0., 1.),
                   linestyle::Symbol = :solid)

    AltPolicy(key, eqcond_fcn, solve, replace_eqcond, forecast_init, color, linestyle)
end

Base.string(a::AltPolicy) = string(a.key)

function altpolicy_replace_eq_entries(key::Symbol)
    if key == :flexible_ait
        return flexible_ait_replace_eq_entries
    elseif key == :zero_rate
        return zero_rate_replace_eq_entries
    elseif key == :smooth_ait_gdp_alt
        smooth_ait_gdp_alt_replace_eq_entries
    elseif key == :smooth_ait_gdp
        smooth_ait_gdp_replace_eq_entries
    elseif key == :ait
        ait_replace_eq_entries
    elseif key == :ngdp
        ngdp_replace_eq_entries
    elseif key == :rw_zero_rate
        rw_zero_rate_replace_eq_entries
    elseif key == :rw
        rw_replace_eq_entries
    else
        error("AltPolicy $(key) has not been added to `altpolicy_replace_eq_entries` as having a method for replacing equilibrium conditions")
    end
end

"""
```
altpolicy_solve(m::AbstractDSGEModel, altpolicy::AltPolicy;
    regime_switching::Bool = false, regimes::Vector{Int} = [1])
```

Solves for the transition equation of `m` under the given altpolicy
"""
function altpolicy_solve(m::AbstractDSGEModel, altpolicy::AltPolicy;
                         regime_switching::Bool = false, regimes::Vector{Int} = [1])
    # Get equilibrium condition matrices

    if length(regimes) == 1
        Γ0, Γ1, C, Ψ, Π  = altpolicy.eqcond(m, regimes[1])
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
            Γ0s[reg], Γ1s[reg], Cs[reg], Ψs[reg], Πs[reg]  = altpolicy.eqcond(m, reg)
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
