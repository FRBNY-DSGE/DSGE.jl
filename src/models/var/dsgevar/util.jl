function draw_stationary_VAR(YYYYD::Matrix{S}, XXYYD::Matrix{S}, XXXXD::Matrix{S},
                             nobsc::Int, nobs::Int, nlags::Int; test::Bool = false,
                             test_Σ_draw::Matrix{S} = Matrix{S}(undef, 0, 0),
                             test_β_draw::Vector{S} = Vector{S}(undef, 0)) where {S<:Real}

    # Set up
    k = 1 + nlags * nobs
    inv_XXXXD = inv(XXXXD)
    β = inv_XXXXD * XXYYD
    inv_Σ_mul_nobsc = inv(YYYYD - XXYYD' * β)
    inv_Σ_mul_nobsc += inv_Σ_mul_nobsc' # force to be positive definite
    inv_Σ_mul_nobsc ./= 2.
    cholmat = cholesky(inv_Σ_mul_nobsc).L

    β = vec(β)
    β_draw = similar(β)
    Σ_draw = similar(inv_Σ_mul_nobsc)
    if test # just do one draw each
        # Draw from marginal posterior of Σ (based on DSGE-VAR)
        z = cholmat * test_Σ_draw
        Σ_draw = inv(z * z')

        # Draw from the conditional posterior of β (based on DSGE-VAR)
        vc       = kron(Σ_draw, inv_XXXXD)
        vc       += vc'
        vc       ./= 2.
        β_draw   = convert(Matrix{S}, reshape(β + cholesky(vc).L * test_β_draw, k, nobs)')
    else
        stationary = false
        while !stationary
            # Draw from marginal posterior of Σ (based on DSGE-VAR)
            z = cholmat * randn(nobs, nobsc - k)
            Σ_draw = inv(z * z')

            # Draw from the conditional posterior of β (based on DSGE-VAR)
            vc       = kron(Σ_draw, inv_XXXXD)
            vc       += vc'
            vc       ./= 2.
            β_draw   = convert(Matrix{S}, reshape(β + cholesky(vc).L * randn(nobs * k), k, nobs)')
            β_to_TTT = vcat(β_draw[:, (1+1):k],
                            hcat(Matrix{S}(I, nobs * (nlags - 1), nobs * (nlags - 1)),
                                 zeros(S, nobs * (nlags - 1), nobs)))
            if maximum(abs.(eigen(β_to_TTT).values)) < 1
                stationary = true
            end
        end
    end

    return β_draw, Σ_draw
end
