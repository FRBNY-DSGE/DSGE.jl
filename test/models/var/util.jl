fp = dirname(@__FILE__)
@testset "Creating lagged data" begin
    # lag
    matdata = matread(joinpath(fp, "../../reference/test_lag.mat"))
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

@testset "Check covariances in data" begin
    matdata = matread(joinpath(fp, "../../reference/test_covariance.mat"))
    lags = Int(matdata["lags"])
    @test DSGE.lag_data(Matrix(matdata["YY"]'), Int(matdata["lags"])) ≈
        matdata["XX"][lags + 1:end, :]
    YYYY, XXYY, XXXX = DSGE.compute_var_covariances(Matrix(matdata["YY"]'), lags)
    @test YYYY ≈ matdata["YY"][lags + 1:end, :]' * matdata["YY"][lags + 1:end, :]
    @test XXYY ≈ matdata["XX"][lags + 1:end, 2:end]' * matdata["YY"][lags + 1:end, :]
    @test XXXX ≈ matdata["XX"][lags + 1:end, 2:end]' * matdata["XX"][lags + 1:end, 2:end]

    # Use intercept now
    YYYY, XXYY, XXXX = DSGE.compute_var_covariances(Matrix(matdata["YY"]'), lags;
                                                    use_intercept = true)
    @test YYYY ≈ matdata["YY"][lags + 1:end, :]' * matdata["YY"][lags + 1:end, :]
    @test XXYY ≈ matdata["XX"][lags + 1:end, 1:end]' * matdata["YY"][lags + 1:end, :]
    @test XXXX ≈ matdata["XX"][lags + 1:end, 1:end]' * matdata["XX"][lags + 1:end, 1:end]
end
