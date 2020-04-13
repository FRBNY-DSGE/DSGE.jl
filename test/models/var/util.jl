fp = dirname(@__FILE__)
@testset "Creating lagged data" begin
    # lag
    matdata = load(joinpath(fp, "../../reference/test_lag.jld2"))
    nonzero = findfirst(vec(sum(matdata["lag_y"], dims = 2)) .!== 0.)
    lags = Int(matdata["lags"])
    @test DSGE.lag(matdata["y"], lags) ≈ matdata["lag_y"][nonzero:end, :]
    @test Matrix(DSGE.lag(Matrix(matdata["y"]'), lags;
              T_by_n = false)') ≈ matdata["lag_y"][nonzero:end, :]
    @test DSGE.lag(matdata["y"], lags;
              drop_obs = lags ÷ 2) ≈
                  matdata["lag_y"][nonzero + (lags ÷ 2):end, :]
    @test Matrix(DSGE.lag(Matrix(matdata["y"]'), lags;
                     drop_obs = lags ÷ 2,
                     T_by_n = false)') ≈ matdata["lag_y"][nonzero + (lags ÷ 2):end, :]

    # lag with padding
    padded = DSGE.lag(matdata["y"], lags, pad = true)
    padded_dropobs1 = DSGE.lag(matdata["y"], lags, pad = true,
                         drop_obs = lags ÷ 2)
    padded_dropobs2 = DSGE.lag(matdata["y"], lags, pad = true,
                          drop_obs = lags + 1)
    padded_dropobs3 = Matrix(DSGE.lag(Matrix(matdata["y"]'), lags, pad = true,
                         drop_obs = lags ÷ 2, T_by_n = false)')
    padded_dropobs4 = Matrix(DSGE.lag(Matrix(matdata["y"]'), lags, pad = true,
                          drop_obs = lags + 1, T_by_n = false)')
    padded[isnan.(padded)] .= 0.
    padded_dropobs1[isnan.(padded_dropobs1)] .= 0.
    padded_dropobs3[isnan.(padded_dropobs3)] .= 0.
    @test padded ≈ matdata["lag_y"]
    @test padded_dropobs1 ≈ matdata["lag_y"][(lags ÷ 2) + 1:end, :]
    @test padded_dropobs2 ≈ matdata["lag_y"][(lags + 2):end, :]
    @test padded_dropobs3 ≈ matdata["lag_y"][(lags ÷ 2) + 1:end, :]
    @test padded_dropobs4 ≈ matdata["lag_y"][(lags + 2):end, :]


    # lag_data
    @test DSGE.lag_data(Matrix(matdata["y"]'), lags) ≈ matdata["x"][lags + 1:end, :]
    @test DSGE.lag_data(Matrix(matdata["y"]'), lags,
                   pad = true)[lags + 1:end, :] ≈ matdata["x"][lags + 1:end, :]
    tmp = DSGE.lag_data(Matrix(matdata["y"]'), lags, pad = true)
    tmp[isnan.(tmp)] .= 0.
    @test tmp ≈ matdata["x"]
    @test DSGE.lag_data(Matrix(matdata["y"]'), lags,
                   use_intercept = false) ≈ matdata["x"][lags + 1:end, 2:end]
    @test DSGE.lag_data(Matrix(matdata["y"]'), lags, pad = true,
                   padding = zeros(lags, size(matdata["x"], 2)))[lags + 1:end, :] ≈
                       matdata["x"][lags + 1:end, :]
end

@testset "Check calculation of population moments" begin
    matdata = load(joinpath(fp, "../../reference/test_var_population_moments.jld2"))
    lags = Int(matdata["lags"])
    @test DSGE.lag_data(Matrix(matdata["YY"]'), Int(matdata["lags"])) ≈
        matdata["XX"][lags + 1:end, :]
    YYYY, XXYY, XXXX = DSGE.compute_var_population_moments(Matrix(matdata["YY"]'), lags)
    @test YYYY ≈ matdata["YY"][lags + 1:end, :]' * matdata["YY"][lags + 1:end, :]
    @test XXYY ≈ matdata["XX"][lags + 1:end, 2:end]' * matdata["YY"][lags + 1:end, :]
    @test XXXX ≈ matdata["XX"][lags + 1:end, 2:end]' * matdata["XX"][lags + 1:end, 2:end]

    # Use intercept now
    YYYY, XXYY, XXXX = DSGE.compute_var_population_moments(Matrix(matdata["YY"]'), lags;
                                                    use_intercept = true)
    @test YYYY ≈ matdata["YY"][lags + 1:end, :]' * matdata["YY"][lags + 1:end, :]
    @test XXYY ≈ matdata["XX"][lags + 1:end, 1:end]' * matdata["YY"][lags + 1:end, :]
    @test XXXX ≈ matdata["XX"][lags + 1:end, 1:end]' * matdata["XX"][lags + 1:end, 1:end]
end

fp = dirname(@__FILE__)
@testset "Draw stationary VAR coefficients and variance-covariance matrix" begin
    matdata = load(joinpath(fp, "../../reference/test_dsgevar_stationary_draws.jld2"))
    β_draw, Σ_draw = DSGE.draw_stationary_VAR(matdata["yyyyd"], matdata["xxyyd"],
                                         matdata["xxxxd"], Int(matdata["nobsc"]),
                                         Int(matdata["nvar"]), Int(matdata["nlags"]);
                                         testing = true, test_Σ_draw_shock = matdata["draw"],
                                         test_β_draw_shock = vec(matdata["beta_draw"]))
    @test matdata["cct_sim"] ≈ convert(Matrix{Float64}, β_draw')
    @test matdata["sig_sim"] ≈ Σ_draw
end

@testset "Draw VECM coefficients and variance-covariance matrix" begin
    matdata = load(joinpath(fp, "../../reference/draw_vecm.jld2"))
    β_draw, Σ_draw = DSGE.draw_VECM(matdata["yyyyc"], matdata["xxyyc"],
                                    matdata["xxxxc"], Int(matdata["Tbar"]),
                                    Int(matdata["n_obs"]), Int(matdata["nlags"]),
                                    Int(matdata["coint"]);
                                    testing = true, test_Σ_draw_shock = matdata["sigma_draw"],
                                    test_β_draw_shock = vec(matdata["beta_draw"]))
    @test maximum(abs.(matdata["beta_ans"] - β_draw)) < 1e-2
    @test @test_matrix_approx_eq matdata["sigma_ans"] Σ_draw
end
