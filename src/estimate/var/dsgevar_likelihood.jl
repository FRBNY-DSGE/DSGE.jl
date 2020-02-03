# data is assumed an nobs x T+lags matrix,
# where the lags indicate we cut off the data for presampling
function dsgevar_likelihood(m::DSGEVAR, data::Matrix{S};
                            λ::S = get_setting(m, :λ), apply_altpolicy::Bool = false,
                            regime_switching::Bool = false,
                            n_regimes::Int = 2, regime::Int = 1,
                            check_system::Bool = true) where {S<:Real}

    lags = n_lags(m)

    # Get covariances
    YYYY, XXYY, XXXX = compute_var_covariances(data, lags)

    yyyyd, xxyyd, xxxxd = compute_system(m; apply_altpolicy = apply_altpolicy,
                                         regime_switching = regime_switching,
                                         regime = regime, n_regimes = n_regimes,
                                         check_system = check_system,
                                         get_covariances = true)

    return dsgevar_likelihood(YYYY, XXYY, XXXX, yyyyd, xxyyd, xxxxd,
                      size(data, 2) - lags, λ, lags, n_observables(m))
end

function dsgevar_likelihood(YYYY::Matrix{S}, XXYY::Matrix{S},
                            XXXX::Matrix{S}, YYYYD::Matrix{S}, XXYYD::Matrix{S},
                            XXXXD::Matrix{S}, T::Int, λ::S,
                            lags::Int, n_obs::Int) where {S<:Real}

    if isinf(λ)
        β = XXXXD \ XXYYD
        Σ = YYYYD - XXYYD' * β
        loglik = -(n_obs * T / 2) * log(2 * π) - (T / 2) * log(det(Σ)) -
            .5 * sum(diag(inv(Σ) * (YYYY - 2 * β' * XXYY + β' * XXXX * β)))
        loglik = real(loglik)
    else
        # Compute weighted covariances
        λT    = λ * T
        YYYYC = YYYY + λT * YYYYD
        XXYYC = XXYY + λT * XXYYD
        XXXXC = XXXX + λT * XXXXD

        # Dimensions
        n_obsc = T + λT
        k     = 1 + lags * n_obs

        # Compute OLS estimates of VAR
        Σ = (YYYYC - XXYYC' * (XXXXC \ XXYYC)) / n_obsc
        Σd = YYYYD - XXYYD' * (XXXXD \ XXYYD)

        # Add terms in normalization constant independent of θ
        loglik = (-n_obs * T / 2) * log(2 * π) + (T * n_obs / 2) * log(2)
        for i = 1:n_obs
            loglik += loggamma((n_obsc - k + 1 - i) / 2) - loggamma((λT - k + 1 - i) / 2)
        end

        # Calculate remainder of likelihood
        # p(Y | θ, λT) = ∫ p(Y | β, Σ) p(β, Σ | θ, λT) d(β, Σ), where
        # - p(Y | β, Σ) is the standard Gaussian likelihood of a VAR unconditional on θ,
        #               i.e. the likelihood ignoring the restrictions a DSGE places on a VAR,
        #                    or the likelihood when doing OLS on just the data.
        # - p(β, Σ | θ, λT) is the prior on (β, Σ) implied by a DSGE-VAR
        #                  with weighting λT.
        loglik += -((n_obsc - k) / 2) * log(det(n_obsc * Σ)) + ((λT - k) / 2) * log(det(λT * Σd)) -
            (n_obs / 2) * log(det(XXXXC) / det(λT * XXXXD))
        loglik = real(loglik)

        # Note that the posterior for DSGE parameters and weight λ is
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
