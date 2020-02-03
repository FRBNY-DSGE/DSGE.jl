# data is assumed an nobs x T+lags matrix,
# where the lags indicate we cut off the data for presampling
function dsgevar_likelihood(m::AbstractDSGEModel, data::Matrix{S};
                            observables::Vector{Symbol} = get_setting(m, :dsgevar_observables),
                            shocks::Vector{Symbol} = collect(keys(m.exogenous_shocks)),
                            lags::Int = get_setting(m, :dsgevar_lags),
                            λ::S = get_setting(m, :dsgevar_λ),
                            regime_switching::Bool = false,
                            n_regimes::Int = 2, regime::Int = 1,
                            MM::Matrix{S} =
                            zeros(length(observables), length(shocks)),
                            include_constant::Bool = false) where {S<:Real}

    YYYY, XXYY, XXXX = compute_var_covariances(data, lags)

    yyyyd, xxyyd, xxxxd = dsgevar(m, observables, shocks, lags;
                                  MM = MM, get_VAR = false)
    # yyyyd, xxyyd, xxxxd = dsgevar(m, observables, shocks, lags;
    #                               regime_switching = regime_switching,
    #                               n_regimes = n_regimes, regime = regime,
    #                               MM = MM, get_VAR = false)
    return dsgevar_likelihood(YYYY, XXYY, XXXX, yyyyd, xxyyd, xxxxd,
                      size(data, 2) - lags, λ, lags, length(observables))
end

function dsgevar_likelihood(YYYY::Matrix{S}, XXYY::Matrix{S},
                            XXXX::Matrix{S}, YYYYD::Matrix{S}, XXYYD::Matrix{S},
                            XXXXD::Matrix{S}, T::Int, λ::S,
                            nlags::Int, nobs::Int) where {S<:Real}

    if isinf(λ)
        β = XXXXD \ XXYYD
        Σ = YYYYD - XXYYD' * β
        loglik = -(nobs * T / 2) * log(2 * π) - (T / 2) * log(det(Σ)) -
            .5 * sum(diag(inv(Σ) * (YYYY - 2 * β' * XXYY + β' * XXXX * β)))
        loglik = real(loglik)
    else
        # Compute weighted covariances
        λT    = λ * T
        YYYYC = YYYY + λT * YYYYD
        XXYYC = XXYY + λT * XXYYD
        XXXXC = XXXX + λT * XXXXD

        # Dimensions
        nobsc = T + λT
        k     = 1 + nlags * nobs

        # Compute OLS estimates of VAR
        Σ = (YYYYC - XXYYC' * inv(XXXXC) * XXYYC) / nobsc
        # Σ = (YYYYC - XXYYC' * (XXXXC \ XXYYC)) / nobsc
        Σd = YYYYD - XXYYD' * inv(XXXXD) * XXYYD
        # Σd = (YYYYD - XXYYD' * (XXXXD \ XXYYD)) / nobsc

        # Add terms in normalization constant independent of θ
        loglik = (-nobs * T / 2) * log(2 * π) + (T * nobs / 2) * log(2)
        for i = 1:nobs
            loglik += loggamma((nobsc - k + 1 - i) / 2) - loggamma((λT - k + 1 - i) / 2)
        end

        # Calculate remainder of likelihood
        # p(Y | θ, λT) = ∫ p(Y | β, Σ) p(β, Σ | θ, λT) d(β, Σ), where
        # - p(Y | β, Σ) is the standard Gaussian likelihood of a VAR unconditional on θ,
        #               i.e. the likelihood ignoring the restrictions a DSGE places on a VAR,
        #                    or the likelihood when doing OLS on just the data.
        # - p(β, Σ | θ, λT) is the prior on (β, Σ) implied by a DSGE-VAR
        #                  with weighting λT.
        loglik += -((nobsc - k) / 2) * log(det(nobsc * Σ)) + ((λT - k) / 2) * log(det(λT * Σd)) -
            (nobs / 2) * log(det(XXXXC) / det(λT * XXXXD))
        loglik = real(loglik)

        # Posterior for DSGE parameters and weight λ is
        # p(θ | Y, λ) ∝ p(Y | θ, λ) p(θ | λ)
        # and joint posterior for VAR parameters is
        # p(β, Σ, θ | Y, λ) = p(β, Σ | θ, Y, λ) p(θ | Y, λ), where
        # p(β, Σ | θ, Y, λ) = p(Y, β, Σ | θ, λ) p(β, Σ | θ, λ).
        # Observe that p(β, Σ | θ, λ) is the prior on (β, Σ) implied by  DSGE-VAR
        #                             with weighting λ, and
        #              p(Y, β, Σ | θ, λ) = p(Y | β, Σ) * (p(β, Σ | θ, λ) * Jeffrey's prior),
        #                                and the right two terms form the Minnesota prior.
    end

    return loglik
end

function impulse_responses(TTT::Matrix{S}, RRR::Matrix{S}, ZZ::Matrix{S},
                           DD::Vector{S}, MM::Matrix{S}, impact::Matrix{S}, horizon::Int;
                           accumulate::Bool = false,
                           cum_inds::Union{Int,UnitRange{Int},Vector{Int}} = 0) where {S<:Real}
    # Get dimensions
    nobs   = size(impact, 1)
    nstate = size(TTT, 1)

    # Compute impulse response to identified impact matrix
    irf_results = Matrix{S}(undef, horizon, nobs^2)
    for i = 1:nobs
        imp    = impact[:, i]
        states = zeros(S, nstate, horizon)
        obs    = zeros(S, nobs, horizon)

        states[:, 1] = RRR * imp
        obs[:, 1]    = ZZ * states[:, 1] + MM * imp
        for t = 2:horizon
            states[:, t] = TTT * states[:, t - 1]
            obs[:, t]    = ZZ * states[:, t] + DD
        end
        if accumulate
            obs[cum_inds, :] = cumsum(obs[cum_inds, :], dims = 2)
        end
        irf_results[:, 1 + (i - 1) * nobs:i * nobs] = obs'
    end

    return irf_results
end
