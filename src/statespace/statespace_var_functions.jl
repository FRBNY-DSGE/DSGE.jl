# This script holds VAR related functions, e.g. DSGE-VAR and DSGE-VECM

"""
```
compute_system(m; apply_altpolicy = false,
            check_system = false, get_system = false,
            get_population_moments = false, use_intercept = false,
            tvis::Bool = false, verbose = :high)
            compute_system(m, data; apply_altpolicy = false,
            check_system = false, get_system = false,
            get_population_moments = false,
            tvis::Bool = false, verbose = :high)
```
Given the current model parameters, compute the DSGE-VAR or DSGE-VECM system
corresponding to model `m`. If a matrix `data` is also passed, then
the VAR is estimated on `data` using the DSGE `m` as a prior
with weight λ.

### Keyword Arguments
* `check_system::Bool`: see `?compute_system` that takes the input `m::AbstractDSGEModel`
   and `system::System`.
* `get_system::Bool`: see Outputs
* `get_population_moments::Bool`: see Outputs
* `use_intercept::Bool`: use an intercept term when computing the OLS estimate of the VAR system.
* `tvis::Bool` indicates whether the state-space system involves time-varying information sets.

### Outputs
* If `get_system = true`:
  Returns the updated `system` whose measurement matrices `ZZ`, `DD`, and `QQ` correspond
  to the VAR or VECM specified by `m`. If `m` is an `AbstractDSGEVECMModel`,
  then the `system` and the vector implied by additional cointegrating relationships
  are returned as a 2-element tuple.
* If `get_population_moments = true`:
  Returns the limit cross product matrices that describe the DSGE implied
    population moments between the observables and their lags. If `data` is
    also passed as an input, then the sample population moments are also returned.
* Otherwise:
Returns `β` and `Σ`, the coefficients and observables covariance matrix of the VAR or VECM.
If `data` is passed in, then `β` and `Σ` are estimated from the data using `m`
as a prior with weight λ. Otherwise, `β` and `Σ` comprise the VECM approximation
of the DSGE `m`.
"""
function compute_system(m::AbstractDSGEVARModel{T}; apply_altpolicy::Bool = false,
                        check_system::Bool = false, get_system::Bool = false,
                        get_population_moments::Bool = false, use_intercept::Bool = false,
                        tvis::Bool = false, verbose::Symbol = :high) where {T <: Real}

    regime_switching = haskey(get_settings(m), :regime_switching) ?
    get_setting(m, :regime_switching) : false
    n_regimes        = regime_switching && haskey(get_settings(m), :n_regimes) ?
    get_setting(m, :n_regimes) : 1

    dsge = get_dsge(m)
    if regime_switching
        error("Regime switching has not been implemented for a DSGEVAR yet.")
        system = compute_system(dsge;
                                verbose = verbose) # This `system` is really a RegimeSwitchingSystem

        systems = Vector{System{T}}(undef, n_regimes)
        for i in 1:n_regimes
            systems[i] = compute_system(dsge, System(system, i); observables = collect(keys(get_observables(m))),
                                        shocks = collect(keys(get_shocks(m))), check_system = check_system)
        end
        system = RegimeSwitchingSystem(systems) # construct a RegimeSwitchingSystem from underlying systems
    else
        system = compute_system(dsge; verbose = verbose)
        system = compute_system(dsge, system; observables = collect(keys(get_observables(m))),
                                shocks = collect(keys(get_shocks(m))), check_system = check_system)
    end

    if get_system
        return system
    elseif regime_switching
        EEs, MMs = measurement_error(m; regime_switching = regime_switching, n_regimes = n_regimes)
        out = get_population_moments ? Vector{Tuple{3, Matrix{T}}}(undef, n_regimes) :
        Vector{Tuple{2, Matrix{T}}}(undef, n_regimes)

        for i in 1:n_regimes
            out[i] = var_approx_state_space(system[i, :TTT], system[i, :RRR], system[i, :QQ],
                                            system[i, :DD], system[i, :ZZ], EEs[i], MMs[i], n_lags(m);
                                            get_population_moments = get_population_moments,
                                            use_intercept = use_intercept)
        end

        return out
    else
        EE, MM = measurement_error(m)

        return var_approx_state_space(system[:TTT], system[:RRR], system[:QQ],
                                      system[:DD], system[:ZZ], EE, MM, n_lags(m);
                                      get_population_moments = get_population_moments,
                                      use_intercept = use_intercept)
    end
end

function compute_system(m::AbstractDSGEVARModel{T}, data::Matrix{T};
                        apply_altpolicy::Bool = false,
                        check_system::Bool = false, get_system::Bool = false,
                        get_population_moments::Bool = false,
                        tvis::Bool = false, verbose::Symbol = :high) where {T<:Real}

    if get_λ(m) == Inf
        # Then we just want the VAR approximation of the DSGE
        return compute_system(m; apply_altpolicy = apply_altpolicy,
                              check_system = check_system, get_system = get_system,
                              get_population_moments = get_population_moments, use_intercept = true,
                              verbose = verbose)
    else
        # Create a system using the method for DSGE-VARs and λ = ∞
        system = compute_system(m; apply_altpolicy = apply_altpolicy,
                                check_system = check_system,
                                get_system = true, use_intercept = true,
                                verbose = verbose)

        if get_system
            return system
        else
            EE, MM = measurement_error(m)

            lags = n_lags(m)
            YYYY, XXYY, XXXX =
            compute_var_population_moments(data, lags; use_intercept = true)
            out = var_approx_state_space(system[:TTT], system[:RRR], system[:QQ],
                                         system[:DD], system[:ZZ], EE, MM, n_lags(m);
                                         get_population_moments = true,
                                         use_intercept = true)

            if get_population_moments
                return out..., YYYY, XXYY, XXXX
            else
                # Compute prior-weighted population moments
                λ = get_λ(m)
                YYYYC = YYYY + λ .* out[1]
                XXYYC = XXYY + λ .* out[2]
                XXXXC = XXXX + λ .* out[3]

                # Draw stationary VAR system
                n_periods = size(data, 2) - lags
                β, Σ =  draw_stationary_VAR(YYYYC, XXYYC, XXXXC,
                                            convert(Int, floor(n_periods + λ * n_periods)),
                                            size(data, 1), lags)

                return β, Σ
            end
        end
    end
end

# Same functions as above but for AbstractDSGEVECMModel types
function compute_system(m::AbstractDSGEVECMModel{T}; apply_altpolicy::Bool = false,
                        check_system::Bool = false, get_system::Bool = false,
                        get_population_moments::Bool = false, use_intercept::Bool = false,
                        tvis::Bool = false, verbose::Symbol = :high) where {T<:Real}

    dsge = get_dsge(m)
    system = compute_system(dsge; verbose = verbose)

    # Use wrapper compute_system for AbstractDSGEVECMModel types
    # (as opposed to compute_system(m::AbstractDSGEModel, system::System; ...))
    system, DD_coint_add = compute_system(m, system; observables = collect(keys(get_observables(m))),
                                          cointegrating = collect(keys(get_cointegrating(m))),
                                          cointegrating_add = collect(keys(get_cointegrating_add(m))),
                                          shocks = collect(keys(get_shocks(m))), check_system = check_system,
                                          get_DD_coint_add = true)

    if get_system
        return system, DD_coint_add
    else
        EE, MM = measurement_error(m)

        return vecm_approx_state_space(system[:TTT], system[:RRR], system[:QQ],
                                       system[:DD], system[:ZZ], EE, MM, n_observables(m),
                                       n_lags(m), n_cointegrating(m), n_cointegrating_add(m),
                                       DD_coint_add;
                                       get_population_moments = get_population_moments,
                                       use_intercept = use_intercept)
    end
end

function compute_system(m::AbstractDSGEVECMModel{T}, data::Matrix{T};
                        apply_altpolicy::Bool = false,
                        check_system::Bool = false, get_system::Bool = false,
                        get_population_moments::Bool = false,
                        tvis::Bool = false, verbose::Symbol = :high) where {T<:Real}

    if get_λ(m) == Inf
        # Then we just want the VECM approximation of the DSGE
        # with no additional cointegration
        return compute_system(m; apply_altpolicy = apply_altpolicy,
                              check_system = check_system, get_system = get_system,
                              get_population_moments = get_population_moments, use_intercept = true,
                              verbose = verbose)
    else
        # Create a system using the method for DSGE-VECMs and λ = ∞
        system, DD_coint_add = compute_system(m; apply_altpolicy = apply_altpolicy,
                                              check_system = check_system,
                                              get_system = true, use_intercept = true,
                                              verbose = verbose)

        if get_system
            return system
        else
            EE, MM = measurement_error(m)

            lags = n_lags(m)
            YYYY, XXYY, XXXX =
            compute_var_population_moments(data, lags; use_intercept = true)
            out = vecm_approx_state_space(system[:TTT], system[:RRR], system[:QQ],
                                          system[:DD], system[:ZZ], EE, MM, size(data, 1),
                                          n_lags(m), n_cointegrating(m),
                                          n_cointegrating_add(m), DD_coint_add;
                                          get_population_moments = true,
                                          use_intercept = true)

            if get_population_moments
                return out..., YYYY, XXYY, XXXX
            else
                # Compute prior-weighted population moments
                λ = get_λ(m)
                YYYYC = YYYY + λ .* out[1]
                XXYYC = XXYY + λ .* out[2]
                XXXXC = XXXX + λ .* out[3]

                # Draw VECM system
                n_periods = size(data, 2) - lags
                β, Σ =  draw_VECM(YYYYC, XXYYC, XXXXC,
                                  convert(Int, n_periods + λ * n_periods),
                                  size(data, 1), lags, n_cointegrating(m))

                return β, Σ
            end
        end
    end
end

"""
```
compute_system(m::PoolModel{T})
```

Given the current model parameters, compute the state-space system
corresponding to the PoolModel model `m`.

Outputs

```
Φ: state transition function
Ψ: likelihood function, given weights on underlying models (the states) and predictive densities
F_ϵ: structural shock distribution
F_u: likelihood function measurement error distribution
F_λ: initial distribution of λ for state transition function
"""
function compute_system(m::PoolModel{T}; tvis::Bool = false, verbose::Symbol = :high) where {T <: Real}
    Φ, F_ϵ, F_λ = transition(m)
    Ψ, F_u = measurement(m)
    return Φ, Ψ, F_ϵ, F_u, F_λ
end

"""
```
compute_system(m::AbstractDSGEModel, system::System;
observables::Vector{Symbol} = collect(keys(m.observables)),
pseudo_observables::Vector{Symbol} = collect(keys(m.pseudo_observables)),
states::Vector{Symbol} = vcat(collect(keys(m.endogenou_states)),
collect(keys(m.endogenous_states_augmented)))
shocks::Vector{Symbol} = collect(keys(m.exogenous_shocks)),
zero_DD = false, zero_DD_pseudo = false)
compute_system(m::AbstractDSGEVECMModel, system::System;
observables::Vector{Symbol} = collect(keys(m.observables)),
pseudo_observables::Vector{Symbol} = collect(keys(m.pseudo_observables)),
cointegrating::Vector{Symbol} = collect(keys(m.cointegrating)),
states::Vector{Symbol} = vcat(collect(keys(m.endogenou_states)),
collect(keys(m.endogenous_states_augmented)))
shocks::Vector{Symbol} = collect(keys(m.exogenous_shocks)),
zero_DD = false, zero_DD_pseudo = false)
```
computes the corresponding transition and measurement
equations specified by the keywords (e.g. `states`, `pseudo_observables`)
using the existing `ZZ`, `ZZ_pseudo`, and `TTT` matrices in `system`.

Note that this function does not update the EE matrix, which is
set to all zeros. To incorporate measurement errors, the user
must specify the EE matrix after applying compute_system.

### Keywords
* `observables`: variables that should be
entered into the new `ZZ` and `DD` matrices.
The `observables` can be both Observables and PseudoObservables,
but they must be an element of the system already.
* `pseudo_observables`: variables that should be
entered into the new `ZZ_pseudo` and `DD_pseudo` matrices.
The `observables` can be both Observables and PseudoObservables,
but they must be an element of the system already.
* `cointegrating`: variables that should be
entered into the new `ZZ` and `DD` matrices as cointegrating relationships.
The `observables` can be both Observables and PseudoObservables,
but they must be an element of the system already.
* `states`: variables that should be
entered into the new `TTT` and `RRR` matrices as states.
They must be existing states.
* `shocks`: variables that should be
entered into the new `RRR` and `QQ` matrices as shocks.
They must be existing exogenous shocks.
"""
function compute_system(m::AbstractDSGEModel{S}, system::System;
                        observables::Vector{Symbol} = collect(keys(m.observables)),
                        pseudo_observables::Vector{Symbol} =
                        collect(keys(m.pseudo_observables)),
                        states::Vector{Symbol} =
                        vcat(collect(keys(m.endogenous_states)),
                             collect(keys(m.endogenous_states_augmented))),
                        shocks::Vector{Symbol} = collect(keys(m.exogenous_shocks)),
                        zero_DD::Bool = false, zero_DD_pseudo::Bool = false,
                        check_system::Bool = true) where {S<:Real}

    # Set up indices
    oid  = m.observables # observables indices dictionary
    pid  = m.pseudo_observables # pseudo observables indices dictionary
    sid  = m.endogenous_states
    said = m.endogenous_states_augmented # states augmented

    # Find shocks to keep
    shock_inds = map(k -> m.exogenous_shocks[k], shocks)
    Qout = system[:QQ][shock_inds, shock_inds]

    # Find states to keep
    if !issubset(states, vcat(collect(keys(sid)), collect(keys(said))))
        false_states = setdiff(states, vcat(collect(keys(sid)), collect(keys(said))))
        error("The following states in keyword `states` do not exist in `system`: " *
              join(string.(false_states), ", "))
    elseif !isempty(setdiff(vcat(collect(keys(sid)), collect(keys(said))), states))
        which_states = Vector{Int}(undef, length(states))
        for i = 1:length(states)
            which_states[i] = haskey(sid, states[i]) ? sid[states[i]] : said[states[i]]
        end
        Tout = system[:TTT][which_states, which_states]
        Rout = system[:RRR][which_states, shock_inds]
        Cout = system[:CCC][which_states]
    else
        which_states = 1:n_states_augmented(m)
        Tout = copy(system[:TTT])
        Rout = system[:RRR][:, shock_inds]
        Cout = copy(system[:CCC])
    end

    # Compute new ZZ and DD matrices if different observables than current system
    if !isempty(symdiff(observables, collect(keys(oid))))
        Zout = zeros(S, length(observables), size(Tout, 1))
        Dout = zeros(S, length(observables))
        for (i, obs) in enumerate(observables)
            Zout[i, :], Dout[i] = if haskey(oid, obs)
                system[:ZZ][oid[obs], which_states], zero_DD ? zero(S) : system[:DD][oid[obs]]
            elseif haskey(pid, obs)
                system[:ZZ_pseudo][pid[obs], which_states], zero_DD ? zero(S) : system[:DD_pseudo][pid[obs]]
            else
                error("Observable/PseudoObservable $obs cannot be found in the DSGE model $m")
            end
        end
    else
        Zout = copy(system[:ZZ])[:, which_states]
        Dout = zero_DD ? zeros(S, size(Zout, 1)) : copy(system[:DD])
    end

    Eout = zeros(S, length(observables), length(observables)) # measurement errors are set to zero

    # Compute new ZZ_pseudo, DD_pseudo if different pseudo_observables than current system
    if !isempty(symdiff(pseudo_observables, collect(keys(pid))))
        Zpseudoout = zeros(S, length(pseudo_observables), size(Tout, 1))
        Dpseudoout = zeros(S, length(pseudo_observables))
        for (i, pseudoobs) in enumerate(pseudo_observables)
            Zpseudoout[i, :], Dpseudoout[i] = if haskey(oid, pseudoobs)
                system[:ZZ][oid[pseudoobs], which_states], zero_DD_pseudo ?
                zero(S) : system[:DD][oid[pseudoobs]]
            elseif haskey(pid, pseudoobs)
                system[:ZZ_pseudo][pid[pseudoobs], which_states], zero_DD_pseudo ?
                zero(S) : system[:DD_pseudo][pid[pseudoobs]]
            else
                error("Observable/PseudoObservable $pseudoobs cannot be found in the DSGE model $m")
            end
        end
    else
        Zpseudoout = copy(system[:ZZ_pseudo])[:, which_states]
        Dpseudoout = zero_DD_pseudo ? zeros(S, size(Zpseudoout, 1)) : copy(system[:DD_pseudo])
    end

    if check_system
        @assert size(Zout, 2) == size(Tout, 1) "Dimension 2 of new ZZ ($(size(Zout,2))) does not match dimension of states ($(size(Tout,1)))."
        @assert size(Qout, 1) == size(Rout, 2) "Dimension 2 of new RRR ($(size(Zout,2))) does not match dimension of shocks ($(size(Qout,1)))."
    end

    # Construct new system
    return System(Transition(Tout, Rout, Cout),
                  Measurement(Zout, Dout, Qout, Eout),
                  PseudoMeasurement(Zpseudoout, Dpseudoout))
end

function compute_system(m::AbstractDSGEVECMModel{S}, system::System;
                        observables::Vector{Symbol} = collect(keys(get_observables(get_dsge(m)))),
                        cointegrating::Vector{Symbol} = Vector{Symbol}(undef, 0),
                        cointegrating_add::Vector{Symbol} = Vector{Symbol}(undef, 0),
                        pseudo_observables::Vector{Symbol} =
                        collect(keys(get_dsge(m).pseudo_observables)),
                        states::Vector{Symbol} =
                        vcat(collect(keys(get_dsge(m).endogenous_states)),
                             collect(keys(get_dsge(m).endogenous_states_augmented))),
                        shocks::Vector{Symbol} = collect(keys(get_dsge(m).exogenous_shocks)),
                        zero_DD::Bool = false, zero_DD_pseudo::Bool = false,
                        get_DD_coint_add::Bool = false,
                        check_system::Bool = true) where {S<:Real}
    # Cointegrating relationships should exist as observables/pseudo_observables already
    # in the underlying DSGE. We assume cointegrating relationships come after normal observables.
    # Default behavior is to recreate the underlying DSGE's state space representation, however.
    sys = compute_system(get_dsge(m), system; observables = vcat(observables, cointegrating),
                         pseudo_observables = pseudo_observables,
                         states = states, shocks = shocks, zero_DD = zero_DD,
                         zero_DD_pseudo = zero_DD_pseudo, check_system = check_system)
    if get_DD_coint_add
        mtype = typeof(m)
        DD_coint_add = if hasmethod(compute_DD_coint_add, (mtype, Vector{Symbol}))
            compute_DD_coint_add(m, cointegrating_add)
        elseif hasmethod(compute_DD_coint_add, (mtype, ))
            compute_DD_coint_add(m)
        else
            compute_DD_coint_add(m, sys, cointegrating_add)
        end
        return sys, DD_coint_add
    else
        return sys
    end
end

"""
```
function compute_DD_coint_add(m::AbstractDSGEVECMModel{S}, system::System,
cointegrating_add::Vector{Symbol}) where {S <: Real}
```
computes `DD_coint_add` for a `DSGEVECM` model. This vector
holds additional cointegrating relationships that do not require
changes to the `ZZ` matrix.

### Note
We recommend overloading this function if there are
cointegrating relationships which a user does not want
to add to the underlying DSGE. The function `compute_system`
first checks for a method `compute_DD_coint_add` that takes
a Type tuple of `(model_type, Vector{Symbol})` and then `(model_type, )`
before calling this method.

This function is generally intended to be internal. As an example of
other such functions, `eqcond` must be user-defined but
is primarily used internally and not directly called by the user in a script.
"""
function compute_DD_coint_add(m::AbstractDSGEVECMModel{S}, system::System,
                              cointegrating_add::Vector{Symbol}) where {S <: Real}
    if !isempty(cointegrating_add)
        oid = get_observables(get_dsge(m))
        pid = get_pseudo_observables(get_dsge(m))
        Dout = zeros(S, length(cointegrating_add))
        for (i, obs) in enumerate(cointegrating_add)
            Dout[i] = if haskey(oid, obs)
                system[:DD][oid[obs]]
            elseif haskey(pid, obs)
                system[:DD_pseudo][pid[obs]]
            else
                error("Observable/PseudoObservable $obs cannot be found in the DSGE model $m")
            end
        end
        return Dout
    else
        @warn "No additional cointegrating relationships specified. Returning an empty vector."
        return Vector{S}(undef, 0)
    end
end

"""
```
var_approx_state_space(TTT, RRR, QQQ, DD, ZZ, EE, MM, p; get_population_moments = false,
use_intercept = false) where {S<:Real}
```
computes the VAR(p) approximation of the linear state space system

```
sₜ = TTT * sₜ₋₁ + RRR * ϵₜ,
yₜ = ZZ * sₜ + DD + uₜ,
```
where the disturbances are assumed to follow
```
ϵₜ ∼ 𝒩 (0, QQ),
uₜ = ηₜ + MM * ϵₜ,
ηₜ ∼ 𝒩 (0, EE).
```
The `MM` matrix implies
```
cov(ϵₜ, uₜ) = QQ * MM'.
```

### Outputs
If `get_population_moments = false`:
* `β`: VAR(p) coefficients
* `Σ`: innovations variance-covariance matrix for the VAR(p) representation
```
yₜ = Xₜβ + μₜ
```
where `Xₜ` appropriately stacks the constants and `p` lags of `yₜ`, and `μₜ ∼ 𝒩 (0, Σ)`.

If `get_population_moments = true`: we return the limit cross product matrices.
* `yyyyd`: 𝔼[y,y]
* `XXXXd`: 𝔼[y,X(lag rr)]
* `XXyyd`: 𝔼[X(lag rr),X(lag ll)]

Using these matrices, the VAR(p) representation is given by
```
β = XXXXd \\ XXyyd
Σ = yyyyd - XXyyd' * β
```

The keyword `use_intercept` specifies whether or not to use an
intercept term in the VAR approximation.
"""
function var_approx_state_space(TTT::AbstractMatrix{S}, RRR::AbstractMatrix{S},
                                QQ::AbstractMatrix{S}, DD::AbstractVector{S},
                                ZZ::AbstractMatrix{S}, EE::AbstractMatrix{S},
                                MM::AbstractMatrix{S}, p::Int;
                                get_population_moments::Bool = false,
                                use_intercept::Bool = false) where {S<:Real}

    n_obs = size(ZZ, 1)

    HH = EE + MM * QQ * MM'
    VV = QQ * MM'

    ## Compute p autocovariances

    ## Initialize Autocovariances
    GAMM0 = zeros(S, n_obs ^ 2, p + 1)
    GA0 =  solve_discrete_lyapunov(TTT, RRR * QQ * RRR')
    Gl   = ZZ * GA0 * ZZ' + ZZ * RRR * VV + (ZZ * RRR * VV)' + HH
    GAMM0[:, 1] = vec(Gl)
    TTl = copy(TTT)
    GA0ZZ = GA0 * ZZ'
    RRRVV = RRR * VV
    for l = 1:p
        Gl = ZZ * TTl * (GA0ZZ + RRRVV) # ZZ * (TTl * GA0Z) * ZZ' + ZZ * (TTl * RRR * VV)
        GAMM0[:, l+1] = vec(Gl)
        TTl = TTl * TTT
    end

    ## Create limit cross product matrices
    yyyyd = zeros(S, n_obs, n_obs)
    if use_intercept
        XXXXd = zeros(S, 1 + p * n_obs, 1 + p * n_obs)
        yyXXd = zeros(S, n_obs, 1 + p * n_obs)

        XXXXd[1, 1] = 1.
        XXXXd[1, 2:1 + p * n_obs] = repeat(DD', 1, p) # same as kron(ones(1, p), DD')
        XXXXd[2:1 + p * n_obs, 1] = repeat(DD, p)     # same as kron(ones(p), DD)
        yyXXd[:, 1] = DD
    else
        XXXXd = zeros(S, p * n_obs, p * n_obs)
        yyXXd = zeros(S, n_obs, p * n_obs)
    end
    DDDD = DD * DD'
    yyyyd = reshape(GAMM0[:, 1], n_obs, n_obs) + DDDD

    shift = use_intercept ? 1 : 0 # for constructing XXXXd, to select the right indices
    for rr = 1:p
        # 𝔼[yy,x(lag rr)]
        yyXXd[:, n_obs * (rr - 1) + 1 + shift:n_obs * rr + shift] =
        reshape(GAMM0[:, rr + 1], n_obs, n_obs) + DDDD

        # 𝔼[x(lag rr),x(lag ll)]
        for ll = rr:p
        yyyydrrll = reshape(GAMM0[:, ll - rr + 1], n_obs, n_obs) + DDDD
        XXXXd[n_obs * (rr - 1) + 1 + shift:n_obs * rr + shift,
              n_obs * (ll - 1) + 1 + shift:n_obs * ll + shift] = yyyydrrll
        XXXXd[n_obs * (ll - 1) + 1 + shift:n_obs * ll + shift,
              n_obs * (rr - 1) + 1 + shift:n_obs * rr + shift] = yyyydrrll'
        end
    end

    XXyyd = convert(Matrix{S}, yyXXd')

    if get_population_moments
        return yyyyd, XXyyd, XXXXd
    else
        β = \(XXXXd, XXyyd)
        Σ = yyyyd - XXyyd' * β
        return β, Σ
    end
end

"""
```
vecm_approx_state_space(TTT, RRR, QQQ, DD, ZZ, EE, MM, n_obs, p, n_coint,
n_coint, n_coint_add, DD_coint_add; get_population_moments = false,
use_intercept = false) where {S<:Real}
```
computes the VECM(p) approximation of the linear state space system

```
sₜ = TTT * sₜ₋₁ + RRR * ϵₜ,
yₜ = ZZ * sₜ + DD + uₜ,
```
where the disturbances are assumed to follow
```
ϵₜ ∼ 𝒩 (0, QQ),
uₜ = ηₜ + MM * ϵₜ,
ηₜ ∼ 𝒩 (0, EE).
```
The `MM` matrix implies
```
cov(ϵₜ, uₜ) = QQ * MM'.
```

### Outputs
If `get_population_moments = false`:
* `β`: VECM(p) coefficients. The first `n_coint + n_coint_add`
coefficients for each observable comprise the error correction terms,
while the following `1 + p * n_obs` terms are the VAR terms.
* `Σ`: innovations variance-covariance matrix for the VECM(p) representation
```
Δyₜ = eₜβₑ + Xₜβᵥ + μₜ
```
where `βₑ` are the coefficients for the error correction terms;
`eₜ` are the error correction terms specifying the cointegrating relationships;
`βᵥ` are the coefficients for the VAR terms; `Xₜ` appropriately stacks
the constants and `p` lags of `Δyₜ`; and `μₜ ∼ 𝒩 (0, Σ)`.

Note that the error correction terms satisfy the mapping
`eₜ' = C * yₜ₋₁`, where `C` is a matrix.

If `get_population_moments = true`: we return the limit cross product matrices.
* `yyyyd`: 𝔼[y,y]
* `XXXXd`: 𝔼[y,X(lag rr)]
* `XXyyd`: 𝔼[X(lag rr),X(lag ll)]

Note that in the rows of `XXyyd` and the rows and columns of `XXXXd`,
the cointegrating relationships are stacked above the constants and
lagged `Δyₜ`.

Using these matrices, the VAR(p) representation is given by
```
β = XXXXd \\ XXyyd
Σ = yyyyd - XXyyd' * β,
```
where `β` has dimensions `n_obs × (n_coint + n_coint_add + 1 + p * n_obs)`,
and `Σ` has dimensions `n_obs × n_obs`.

The keyword `use_intercept` specifies whether or not to use an
intercept term in the VECM approximation.
"""
function vecm_approx_state_space(TTT::AbstractMatrix{S}, RRR::AbstractMatrix{S},
                                 QQ::AbstractMatrix{S}, DD::AbstractVector{S},
                                 ZZ::AbstractMatrix{S}, EE::AbstractMatrix{S},
                                 MM::AbstractMatrix{S}, n_obs::Int, p::Int,
                                 n_coint::Int, n_coint_add::Int = 0,
                                 DD_coint_add::AbstractVector{S} = Vector{S}(undef, 0);
                                 get_population_moments::Bool = false,
                                 use_intercept::Bool = false,
                                 test_GA0::AbstractMatrix{S} =
                                 Matrix{S}(undef, 0, 0)) where {S <: Real}

    # n_obs is number of observables, n_coint is number of cointegrating relationships
    n_coint_all = n_coint + n_coint_add

    # Create variance-covariance matrices w/measurement error included
    HH = EE + MM * QQ * MM'
    VV = QQ * MM'

    ## Compute p autocovariances

    ## Initialize Autocovariances
    GAMM0 = zeros(S, (n_obs + n_coint) ^ 2, p + 1)
    GA0 = isempty(test_GA0) ? solve_discrete_lyapunov(TTT, RRR * QQ * RRR') : test_GA0 # Matlab code uses a different numerical procedure -> some difference in this matrix
    Gl   = ZZ * GA0 * ZZ' + ZZ * RRR * VV + (ZZ * RRR * VV)' + HH
    GAMM0[:, 1] = vec(Gl)
    TTl = copy(TTT)
    GA0ZZ = GA0 * ZZ'
    RRRVV = RRR * VV
    for l = 1:p
        Gl = ZZ * TTl * (GA0ZZ + RRRVV) # ZZ * (TTl * GA0) * ZZ' + ZZ * (TTl * RRR * VV)
        GAMM0[:, l+1] = vec(Gl)
        TTl = TTl * TTT
    end

    ## Create limit cross product matrices
    DDDD = DD * DD'
    yyyyd_coint0 = reshape(GAMM0[:, 1], n_obs + n_coint, n_obs + n_coint) + DDDD
    yyyyd_coint1 = reshape(GAMM0[:, 2], n_obs + n_coint, n_obs + n_coint) + DDDD
    yyyyd = yyyyd_coint0[1:n_obs, 1:n_obs]

    # n_coint_add are treated as the first set of variables in XX
    # n_coint    are treated as the second set of variables in XX
    # composition: n_coint_add - n_coint - constant - lags
    if use_intercept
        XXXXd = zeros(S, 1 + p * n_obs + n_coint_all, 1 + p * n_obs + n_coint_all)
        yyXXd = zeros(S, n_obs, 1 + p * n_obs + n_coint_all)

        # 𝔼[x(n_coint), x(n_coint)]
        XXXXd[n_coint_add + 1:n_coint_all, n_coint_add + 1:n_coint_all] = yyyyd_coint0[n_obs + 1:n_obs + n_coint, n_obs + 1:n_obs + n_coint]

        # 𝔼[x(const), x(const)]
        XXXXd[n_coint_all + 1, n_coint_all + 1] = 1.

        # 𝔼[x(n_coint), x(const)]
        XXXXd[n_coint_add + 1:n_coint_all, n_coint_all + 1] = DD[n_obs + 1:n_obs + n_coint]
        XXXXd[n_coint_all + 1, n_coint_add + 1:n_coint_all] = DD[n_obs + 1:n_obs + n_coint]

        # 𝔼[x(const), x(lags)]
        XXXXd[n_coint_all + 1, n_coint_all + 2:n_coint_all + 1 + p * n_obs] = repeat(DD[1:n_obs]', 1, p)
        XXXXd[n_coint_all + 2:n_coint_all + 1 + p * n_obs, n_coint_all + 1] =
        XXXXd[n_coint_all + 1, n_coint_all + 2:n_coint_all + 1 + p * n_obs]' # same as kron(ones(p), DD[1:n_obs]) but avoids calculation

        # 𝔼[yy, x(n_coint)]
        yyXXd[:, n_coint_add + 1:n_coint_all] = yyyyd_coint1[1:n_obs, n_obs + 1:n_obs + n_coint]

        # 𝔼[yy, x(const)]
        yyXXd[:, n_coint_all + 1] = DD[1:n_obs]
    else
        XXXXd = zeros(S, p * n_obs + n_coint_all, p * n_obs + n_coint_all)
        yyXXd = zeros(S, n_obs, p * n_obs + n_coint_all)

        # 𝔼[x(n_coint), x(n_coint)]
        XXXXd[n_coint_add + 1:n_coint_all, n_coint_add + 1:n_coint_all] =
        yyyyd_coint0[n_obs + 1:n_obs + n_coint, n_obs + 1:n_obs + n_coint]

        # 𝔼[yy, x(n_coint)]
        yyXXd[:, n_coint_add + 1:n_coint_all] = yyyyd_coint1[1:n_obs, n_obs + 1:n_obs + n_coint]
    end

    if n_coint_add > 0
        DD_coint_add_div2 = DD_coint_add ./ 2
        if use_intercept
            # 𝔼[yy, x(n_coint_add)]
            yyXXd[:, 1:n_coint_add] = DD[1:n_obs] * DD_coint_add_div2

            # 𝔼[x(n_coint_add), x(n_coint_add)]
            XXXXd[1:n_coint_add, 1:n_coint_add] = DD_coint_add * DD_coint_add' ./ 3

            # 𝔼[x(n_coint_add), x(n_coint)]
            XXXXd[1:n_coint_add, 1 + n_coint_add:n_coint_all] =
            DD_coint_add_div2 * DD[n_obs + 1: n_obs + n_coint]'
            XXXXd[1 + n_coint_add:n_coint_all, 1:n_coint_add] =
            XXXXd[1:n_coint_add, 1 + n_coint_add:n_coint_all]' # transpose of the previous line

            # 𝔼[x(n_coint_add), x(const)]
            XXXXd[1:n_coint_add, n_coint_all + 1] = DD_coint_add_div2
            XXXXd[n_coint_all + 1, 1:n_coint_add] = DD_coint_add_div2'
        else
            # 𝔼[yy, x(n_coint_add)]
            yyXXd[:, 1:n_coint_add] = DD[1:n_obs] * DD_coint_add_div2

            # 𝔼[x(n_coint_add), x(n_coint_add)]
            XXXXd[1:n_coint_add, 1:n_coint_add] = DD_coint_add * DD_coint_add' ./ 3

            # 𝔼[x(n_coint_add), x(n_coint)]
            XXXXd[1:n_coint_add, 1 + n_coint_add:n_coint_all] =
            DD_coint_add_div2 * DD[n_obs + 1: n_obs + n_coint]'
            XXXXd[1 + n_coint_add:n_coint_all, 1:n_coint_add] =
            XXXXd[1:n_coint_add, 1 + n_coint_add:n_coint_all]' # transpose of the previous line
        end
    end

    shift = use_intercept ? 1 : 0 # for constructing XXXXd, to select the right indices
    for rr = 1:p
        # 𝔼[yy, x(lag rr)]
        yyyyd_cointrr = reshape(GAMM0[:, rr + 1], n_obs + n_coint, n_obs + n_coint) + DDDD
        yyXXd[:, n_coint_all + 1 + n_obs * (rr - 1) + shift:n_coint_all + n_obs * rr + shift] =
        yyyyd_cointrr[1:n_obs, 1:n_obs]

        if n_coint_add > 0
            # 𝔼[x(n_coint_add), x(lag rr)]
            XXXXd[1:n_coint_add, n_coint_all + 1 + n_obs * (rr - 1) + shift:n_coint_all + n_obs * rr + shift] =
            DD_coint_add_div2 * DD[1:n_obs]'
            XXXXd[n_coint_all + 1 + n_obs * (rr - 1) + shift:n_coint_all + n_obs * rr + shift, 1:n_coint_add] =
            XXXXd[1:n_coint_add, n_coint_all + 1 + n_obs * (rr - 1) + shift:n_coint_all + n_obs * rr + shift]'
        end

        # 𝔼[x(n_coint), x(lag rr)]
        yyyyd_cointrr1 = reshape(GAMM0[:, rr], n_obs + n_coint, n_obs + n_coint) + DDDD
        XXXXd[n_coint_add + 1:n_coint_all, n_coint_all + 1 + n_obs * (rr - 1) + shift:n_coint_all + n_obs * rr + shift] =
        yyyyd_cointrr1[n_obs + 1:n_obs + n_coint, 1:n_obs]
        XXXXd[n_coint_all + 1 + n_obs * (rr - 1) + shift:n_coint_all + n_obs * rr + shift, n_coint_add + 1:n_coint_all] =
        yyyyd_cointrr1[n_obs + 1:n_obs + n_coint, 1:n_obs]'

        # 𝔼[x(lag rr), x(lag ll)]
        for ll = rr:p
            yyyyd_cointrrll = reshape(GAMM0[:, ll - rr + 1], n_obs + n_coint, n_obs + n_coint) + DDDD
            XXXXd[n_coint_all + 1 + n_obs * (rr - 1) + shift:n_coint_all + n_obs * rr + shift,
                  n_coint_all + 1 + n_obs * (ll - 1) + shift:n_coint_all + n_obs * ll + shift] = yyyyd_cointrrll[1:n_obs, 1:n_obs]
            XXXXd[n_coint_all + 1 + n_obs * (ll - 1) + shift:n_coint_all + n_obs * ll + shift,
                  n_coint_all + 1 + n_obs * (rr - 1) + shift:n_coint_all + n_obs * rr + shift] = yyyyd_cointrrll[1:n_obs, 1:n_obs]'
        end
    end

    XXyyd = convert(Matrix{S}, yyXXd')

    if get_population_moments
        return yyyyd, XXyyd, XXXXd
    else
        β = XXXXd \ XXyyd
        Σ = yyyyd - XXyyd' * β
        return β, Σ
    end
end
