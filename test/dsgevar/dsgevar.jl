using DSGE, MAT, Test, CSV

@testset "Likelihood function of DSGE-VAR p(Y | θ, λ)" begin
    matdata = matread("../reference/test_dsgevar_likelihood.mat")
    @test matdata["lnpy0"] ≈ DSGE.dsgevar_likelihood(matdata["YYYY"],
                                                matdata["XXYY"], matdata["XXXX"],
                                                matdata["yyyyd"], matdata["xxyyd"],
                                                matdata["xxxxd"], Int(matdata["nobs"]),
                                                matdata["lambda"] / matdata["nobs"],
                                                Int(matdata["nlags"]),
                                                Int(matdata["nvar"]))
    @test matdata["lnpyinf"] ≈ DSGE.dsgevar_likelihood(matdata["YYYY"],
                                                  matdata["XXYY"], matdata["XXXX"],
                                                  matdata["yyyyd"], matdata["xxyyd"],
                                                  matdata["xxxxd"], Int(matdata["nobs"]),
                                                  Inf, Int(matdata["nlags"]),
                                                  Int(matdata["nvar"]))


    m = Model1002("ss10")
    obs_i = [m.observables[:obs_nominalrate], m.observables[:obs_gdp], m.observables[:obs_gdpdeflator]]
    observables = [:obs_nominalrate, :obs_gdp, :π_t]
    dsge_data = df_to_matrix(m, CSV.read("../reference/test_dsgevar_likelihood_dsge_data.csv"))[obs_i, :]
    @test DSGE.dsgevar_likelihood(m, dsge_data; observables = observables, lags = 4, λ = .5 * size(dsge_data, 2)) !== 0.
end

@testset "Draw stationary VAR coefficients and covariance matrix" begin
    matdata = matread("../reference/test_dsgevar_stationary_draws.mat")
    β_draw, Σ_draw = DSGE.draw_stationary_VAR(matdata["yyyyd"], matdata["xxyyd"],
                                         matdata["xxxxd"], Int(matdata["nobsc"]),
                                         Int(matdata["nvar"]), Int(matdata["nlags"]);
                                         test = true, test_Σ_draw = matdata["draw"],
                                         test_β_draw = vec(matdata["beta_draw"]))
    @test matdata["cct_sim"] ≈ β_draw
    @test matdata["sig_sim"] ≈ Σ_draw
end

@testset "Impulse responses to identified impact matrix" begin
    matdata = matread("../reference/test_irfdsge.mat")
    irfout = impulse_responses(matdata["TTT"], matdata["RRR"], matdata["zz"],
                               zeros(size(matdata["zz"], 1)), matdata["mmm"],
                               matdata["impact"], 1; accumulate = true,
                               cum_inds = 1)
    @test matdata["aairf"] ≈ irfout
end

@testset "Rotation of impulse responses to identified impact matrix" begin
    matdata1 = matread("../reference/test_irfdsge.mat")
    matdata2 = matread("../reference/test_rotations.mat")
    ŷ = DSGE.impulse_responses_rotation(matdata1["TTT"], matdata1["RRR"], matdata1["zz"],
                                        zeros(size(matdata1["zz"], 1)), matdata1["mmm"],
                                        matdata1["impact"], Int(matdata2["k"]), matdata2["cct_sim"],
                                        matdata2["sig_sim"], matdata2["XXpred"],
                                        Int(matdata2["qahead"]); accumulate = true,
                               cum_inds = 1, test_shocks = matdata2["Shocks"])

    @test ŷ ≈ matdata2["yypred"]
end

@testset "Creating lagged data" begin
    # lag
    matdata = matread("../reference/test_lag.mat")
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
                   include_constant = false) ≈ matdata["x"][lags + 1:end, 2:end]
    @test DSGE.lag_data(Matrix(matdata["y"]'), lags, pad = true,
                   padding = zeros(lags, size(matdata["x"], 2)))[lags + 1:end, :] ≈
                       matdata["x"][lags + 1:end, :]
end

@testset "Check covariances in data" begin
    matdata = matread("../reference/test_covariance.mat")
    lags = Int(matdata["lags"])
    @test DSGE.lag_data(Matrix(matdata["YY"]'), Int(matdata["lags"])) ≈
        matdata["XX"][lags + 1:end, :]
    YYYY, XXYY, XXXX = DSGE.compute_var_covariances(Matrix(matdata["YY"]'), lags)
    @test YYYY ≈ matdata["YY"][lags + 1:end, :]' * matdata["YY"][lags + 1:end, :]
    @test XXYY ≈ matdata["XX"][lags + 1:end, 2:end]' * matdata["YY"][lags + 1:end, :]
    @test XXXX ≈ matdata["XX"][lags + 1:end, 2:end]' * matdata["XX"][lags + 1:end, 2:end]
end


nothing
