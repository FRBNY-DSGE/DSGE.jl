function lag(X::Matrix{S}, lags::Int; pad::Bool = false, T_by_n::Bool = true,
             drop_obs::Int = 0) where {S<:Real}
    if pad
        keep_lags = max(lags - drop_obs, 0)
        start = 1 + max(drop_obs - lags, 0)
        if T_by_n
            out = fill!(Matrix{S}(undef, size(X, 1) - drop_obs, size(X, 2)), NaN)
            out[1 + keep_lags:end, :] = X[start:size(X, 1) - lags, :]
        else
            out = fill!(Matrix{S}(undef, size(X, 1), size(X, 2) - drop_obs), NaN)
            out[:, 1 + keep_lags:end] = X[:, start:size(X, 2) - lags]
        end
    else
        out = T_by_n ? X[1 + drop_obs:size(X, 1) - lags, :] :
            X[:, 1 + drop_obs:size(X, 2) - lags]
    end
    return out
end

function lag_data(data::Matrix{S}, lags::Int; use_intercept::Bool = true,
                  pad::Bool = false, padding::Matrix{S} = Matrix{S}(undef, 0, 0)) where {S<:Real}
    # Assumes data is nobs x T
    nobs, T = size(data)
    data = Matrix(data')

    # Construct XX matrix of covariates
    add_constant = use_intercept ? 1 : 0
    XX = fill!(Matrix{S}(undef, pad ? T : T - lags, lags * nobs + add_constant), NaN)
    if use_intercept
        XX[:, 1] .= one(S) # XX is T x n_regressors
    end

    for i = 1:lags
        XX[:, add_constant + (i - 1) * nobs + 1:add_constant + i * nobs] =
            lag(data, i; pad = pad, T_by_n = true, drop_obs = pad ? 0 : lags - i)
    end

    if pad && !isempty(padding)
        XX[1:lags, :] = padding
    end

    return XX
end

function compute_var_population_moments(data::Matrix{S}, lags::Int;
                                 use_intercept::Bool = false) where {S<:Real}
    # Compute population moments of sample data
    YY = convert(Matrix{S}, data[:, 1 + lags:end]')
    XX = lag_data(data, lags; use_intercept = use_intercept) # Construct XX matrix of covariates
    YYYY = YY' * YY
    XXYY = XX' * YY
    XXXX = XX' * XX

    return YYYY, XXYY, XXXX
end


# function draw_stationary_VAR(YYYYD::Matrix{S}, XXYYD::Matrix{S}, XXXXD::Matrix{S},
#                              nobsc::Int, nobs::Int, nlags::Int; test::Bool = false,
#                              test_Σ_draw::Matrix{S} = Matrix{S}(undef, 0, 0),
#                              test_β_draw::Vector{S} = Vector{S}(undef, 0)) where {S<:Real}

#     # Set up
#     k = 1 + nlags * nobs
#     inv_XXXXD = inv(XXXXD)
#     β = inv_XXXXD * XXYYD
#     inv_Σ_mul_nobsc = inv(YYYYD - XXYYD' * β)
#     inv_Σ_mul_nobsc += inv_Σ_mul_nobsc' # force to be positive definite
#     inv_Σ_mul_nobsc ./= 2.
#     cholmat = cholesky(inv_Σ_mul_nobsc).L

#     β = vec(β)
#     β_draw = similar(β)
#     Σ_draw = similar(inv_Σ_mul_nobsc)
#     if test # just do one draw each
#         # Draw from marginal posterior of Σ (based on DSGE-VAR)
#         z = cholmat * test_Σ_draw
#         Σ_draw = inv(z * z')

#         # Draw from the conditional posterior of β (based on DSGE-VAR)
#         vc       = kron(Σ_draw, inv_XXXXD)
#         vc       += vc'
#         vc       ./= 2.
#         β_draw   = convert(Matrix{S}, reshape(β + cholesky(vc).L * test_β_draw, k, nobs)')
#     else
#         stationary = false
#         while !stationary
#             # Draw from marginal posterior of Σ (based on DSGE-VAR)
#             z = cholmat * randn(nobs, nobsc - k)
#             Σ_draw = inv(z * z')

#             # Draw from the conditional posterior of β (based on DSGE-VAR)
#             vc       = kron(Σ_draw, inv_XXXXD)
#             vc       += vc'
#             vc       ./= 2.
#             β_draw   = convert(Matrix{S}, reshape(β + cholesky(vc).L * randn(nobs * k), k, nobs)')
#             β_to_TTT = vcat(β_draw[:, (1+1):k],
#                             hcat(Matrix{S}(I, nobs * (nlags - 1), nobs * (nlags - 1)),
#                                  zeros(S, nobs * (nlags - 1), nobs)))
#             if maximum(abs.(eigen(β_to_TTT).values)) < 1
#                 stationary = true
#             end
#         end
#     end

#     return β_draw, Σ_draw
# end
