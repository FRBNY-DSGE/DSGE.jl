"""
```
dsgevecm_likelihood(m::AbstractDSGEVECMModel{S}, data::Matrix{S};
    apply_altpolicy::Bool = false) where {S<:Real}
```
evaluates the likelihood of a DSGE-VECM. It is assumed that
`get_λ(m)` retrieves the prior weight `λ` on the DSGE.


The matrix `data` is assumed an nobs x T+lags matrix,
where the lags indicate we cut off the data for presampling.

In the future, this function may be combined with `dsgevar_likelihood`
because the two functions are almost exactly the same. We currently have
them implemented separately so that the naming convention does not
create any confusion. One can think of DSGE-VECMs as DSGE-VARs
with an additional regressor that corrects for errors to account for
cointegration, hence it would not be "wrong" per se to add
a couple `if-else` conditions inside `dsgevar_likelihood` and let this
function cover both DSGE-VARs and DSGE-VECMS. But the current decision
to have separate functions makes it clear that the `dsgevar_likelihood`
specifically operates on DSGE-VARs without an error correction term, and
`dsgevecm_likelihood` operates on DSGE-VECMs specifically.
"""
function dsgevecm_likelihood(m::AbstractDSGEVECMModel{S}, data::Matrix{S};
                            apply_altpolicy::Bool = false,
                            use_intercept::Bool = true) where {S<:Real}

    # Set up. Default to non-regime switching if these settings don't exist.
    lags = n_lags(m)
    n_regimes        = haskey(get_settings(m), :n_regimes) ?
        get_setting(m, :n_regimes) : 1
    regime_switching = haskey(get_settings(m), :regime_switching) && n_regimes > 1 ?
        get_setting(m, :regime_switching) : false # if n_regimes == 1, then no switching needed

    # Get population moments
    YYYY, XXYY, XXXX = compute_var_population_moments(data, lags; use_intercept = use_intercept)

    if regime_switching
        error("Regime switching has not been implemented yet.")
        yyyyd = OrderedDict{Int, Matrix{S}}()
        xxyyd = OrderedDict{Int, Matrix{S}}()
        xxxxd = OrderedDict{Int, Matrix{S}}()

        for reg = 1:n_regimes
            yyyyd[reg], xxyyd[reg], xxxxd[reg] =
                compute_system(m; apply_altpolicy = apply_altpolicy,
                               regime_switching = regime_switching,
                               regime = reg, n_regimes = n_regimes,
                               get_population_moments = true,
                               use_intercept = use_intercept)
        end
    else
        yyyyd, xxyyd, xxxxd = compute_system(m; apply_altpolicy = apply_altpolicy,
                                             get_population_moments = true,
                                             use_intercept = use_intercept)
    end

    return dsgevecm_likelihood(YYYY, XXYY, XXXX, yyyyd, xxyyd, xxxxd,
                      size(data, 2) - lags, get_λ(m); n_cointadd = n_cointegrating_add(m))
end

"""
```
dsgevecm_likelihood(YYYY::Matrix{S}, XXYY::Matrix{S},
    XXXX::Matrix{S}, YYYYD::Matrix{S}, XXYYD::Matrix{S},
    XXXXD::Matrix{S}, T::Int, λ::S,
    lags::Int, n_obs::Int,
    coint_inds::Union{Vector{Int}, UnitRange{Int}}) where {S<:Real}
```
evaluates the likelihood of a DSGE-VECM given the
population moments of the raw data (`YYYY`, `XXYY`, `XXXX`)
and the population moments implied by the DSGE prior (`YYYYD`, `XXYYD`, `XXXXD`).

### Other Inputs
* `T`: number of time periods in data
* `λ`: weight placed on the DSGE prior
* `lags`: number of lags in the VECM
* `n_obs`: number of observables (e.g. number of time series).
* `coint_inds`: indices of the cointegrated observables
"""
function dsgevecm_likelihood(YYYY::Matrix{S}, XXYY::Matrix{S},
                             XXXX::Matrix{S}, YYYYD::Matrix{S}, XXYYD::Matrix{S},
                             XXXXD::Matrix{S}, T::Int, λ::S;
                             n_cointadd::Int = 0) where {S <: Real}

    if isinf(λ)
        if n_cointadd > 0
            inds = 1 + n_cointadd:size(XXXX, 1)
            XXXXD = XXXXD[inds, inds]
            XXYYD = XXYYD[inds, :]
            XXXX  = XXXXD[inds, inds]
            XXYY  = XXYY[inds, :]
        end
        β = XXXXD \ XXYYD
        Σ = YYYYD - XXYYD' * β
        n_obs  = size(YYYY, 1)
        loglik = -(n_obs * T / 2) * log(2 * π) - (T / 2) * log(det(Σ)) -
            .5 * sum(diag(Σ \ (YYYY - β' * XXYY - XXYY' * β + β' * XXXX * β)))
        return real(loglik)
    else
        # Compute DSGEVECM components for weighted population moments
        if n_cointadd > 0
            Λ⁻¹ = Diagonal{S}(I, size(XXXX, 1))    # these lines are the only ones which distinguish
            Λ⁻¹[1:n_cointadd, 1:n_cointadd] ./= λT # this likelihood from dsgevar_likelihood
            XXYY = Λ⁻¹ * XXYY
            XXXX = Λ⁻¹ * XXXX * Λ⁻¹
        end

        return dsgevar_likelihood(YYYY, XXYY, XXXX, YYYYD, XXYYD, XXXXD, T, λ)
    end
end
