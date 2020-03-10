fp = dirname(@__FILE__)
@testset "DSGE impulse responses to a pre-specified impact matrix" begin
    matdata = load(joinpath(fp, "../../../reference/test_irfdsge.jld2"))
    irfout = impulse_responses(matdata["TTT"], matdata["RRR"], matdata["zz"],
                               zeros(size(matdata["zz"], 1)), matdata["mmm"],
                               matdata["impact"], 1; accumulate = true, cum_inds = 1)
    nobs = size(matdata["zz"], 1)
    for i = 1:nobs
        @test vec(matdata["aairf"][:, nobs * (i - 1) + 1:nobs * i]) ≈
            vec(irfout[:, :, i])
    end

    irfout_accum1 = impulse_responses(matdata["TTT"], matdata["RRR"], matdata["zz"],
                                      zeros(size(matdata["zz"], 1)), matdata["mmm"],
                                      matdata["impact"], 2; accumulate = true, cum_inds = 1)
    irfout_nocum1 = impulse_responses(matdata["TTT"], matdata["RRR"], matdata["zz"],
                                      zeros(size(matdata["zz"], 1)), matdata["mmm"],
                                      matdata["impact"], 2; accumulate = false, cum_inds = 1)
    irfout_accum2 = impulse_responses(matdata["TTT"], matdata["RRR"], matdata["zz"],
                                      zeros(size(matdata["zz"], 1)), matdata["mmm"],
                                      matdata["impact"], 2; accumulate = true,
                                      cum_inds = 1:2)
    irfout_nocum2 = impulse_responses(matdata["TTT"], matdata["RRR"], matdata["zz"],
                                      zeros(size(matdata["zz"], 1)), matdata["mmm"],
                                      matdata["impact"], 2; accumulate = false,
                                      cum_inds = 1:2)
    irfout_accum3 = impulse_responses(matdata["TTT"], matdata["RRR"], matdata["zz"],
                                      zeros(size(matdata["zz"], 1)), matdata["mmm"],
                                      matdata["impact"], 2; accumulate = true,
                                      cum_inds = [1,3])
    irfout_nocum3 = impulse_responses(matdata["TTT"], matdata["RRR"], matdata["zz"],
                                      zeros(size(matdata["zz"], 1)), matdata["mmm"],
                                      matdata["impact"], 2; accumulate = false,
                                      cum_inds = [1,3])
    for i = 1:size(irfout_accum1, 3)
        @test cumsum(irfout_nocum1[1, :, i]) ≈ irfout_accum1[1, :, i]
        @test irfout_nocum1[2:end, :, i] ≈ irfout_accum1[2:end, :, i]
        @test cumsum(irfout_nocum2[1:2, :, i]; dims = 2) ≈ irfout_accum2[1:2, :, i]
        @test irfout_nocum2[end, :, i] ≈ irfout_accum2[end, :, i]
        @test cumsum(irfout_nocum3[[1,3], :, i]; dims = 2) ≈ irfout_accum3[[1,3], :, i]
        @test irfout_nocum3[2, :, i] ≈ irfout_accum3[2, :, i]
    end
end

@testset "Impulse responses to structural shocks identified by a DSGE-VAR" begin
    matdata1 = load(joinpath(fp, "../../../reference/test_irfdsge.jld2"))
    matdata2 = load(joinpath(fp, "../../../reference/test_rotations.jld2"))
    ŷ = DSGE.impulse_responses(matdata1["TTT"], matdata1["RRR"], matdata1["zz"],
                               zeros(size(matdata1["zz"], 1)), matdata1["mmm"],
                               matdata1["impact"].^2, Int(matdata2["k"]),
                               convert(Matrix{Float64}, matdata2["cct_sim"]'),
                               matdata2["sig_sim"], vec(matdata2["XXpred"]),
                               Int(matdata2["qahead"]); test_shocks =
                               convert(Matrix{Float64}, matdata2["Shocks"]'))

    @test ŷ ≈ matdata2["yypred"]'

    Random.seed!(1793) # drawing shocks, so this shouldn't be the same
    ŷ1 = DSGE.impulse_responses(matdata1["TTT"], matdata1["RRR"], matdata1["zz"],
                                zeros(size(matdata1["zz"], 1)), matdata1["mmm"],
                                matdata1["impact"].^2, Int(matdata2["k"]),
                                convert(Matrix{Float64}, matdata2["cct_sim"]'),
                                matdata2["sig_sim"], vec(matdata2["XXpred"]),
                                Int(matdata2["qahead"]), draw_shocks = true)
    Random.seed!(1793) # Flipping shocks here
    ŷ2 = DSGE.impulse_responses(matdata1["TTT"], matdata1["RRR"], matdata1["zz"],
                                zeros(size(matdata1["zz"], 1)), matdata1["mmm"],
                                matdata1["impact"].^2, Int(matdata2["k"]),
                                convert(Matrix{Float64}, matdata2["cct_sim"]'),
                                matdata2["sig_sim"], vec(matdata2["XXpred"]),
                                Int(matdata2["qahead"]), draw_shocks = true)
    @test !(ŷ1 ≈ matdata2["yypred"]')
    @test ŷ1 ≈ ŷ2
end

@testset "Impulse responses of a VAR using a DSGE as a prior" begin
    jlddata = load(joinpath(fp, "../../../reference/test_dsgevar_lambda_irfs.jld2"))
    m = Model1002("ss10", custom_settings =
                  Dict{Symbol,Setting}(:add_laborshare_measurement =>
                                       Setting(:add_laborshare_measurement, true)))
    m <= Setting(:impulse_response_horizons, 10)
    dsgevar = DSGEVAR(m, collect(keys(m.exogenous_shocks)), "ss11")
    DSGE.update!(dsgevar, λ = 1.)
    DSGE.update!(dsgevar, jlddata["modal_param"])

    Random.seed!(1793)
    out = impulse_responses(dsgevar, jlddata["data"], :cholesky, 1; use_intercept = true)
    out_lr = impulse_responses(dsgevar, jlddata["data"], :choleskyLR, 1; use_intercept = true)
    out_maxbc = impulse_responses(dsgevar, jlddata["data"], :maxBC, 1; use_intercept = true)

    Random.seed!(1793)
    _ = impulse_responses(dsgevar, jlddata["data"], :cholesky, 1; use_intercept = true)
    out_lr2 = impulse_responses(dsgevar, jlddata["data"], :cholesky_long_run, 1;
                                use_intercept = true)
    out_maxbc2 = impulse_responses(dsgevar, jlddata["data"], :maximum_business_cycle_variance,
                                   1; use_intercept = true)

    Random.seed!(1793)
    out_flip = impulse_responses(dsgevar, jlddata["data"], :cholesky, 1; use_intercept = true,
                                 flip_shocks = true,)
    out_lr_flip = impulse_responses(dsgevar, jlddata["data"], :choleskyLR, 1; use_intercept = true,
                                    flip_shocks = true)
    out_maxbc_flip = impulse_responses(dsgevar, jlddata["data"], :maxBC, 1; use_intercept = true,
                                       flip_shocks = true)

    Random.seed!(1793)
    out_h = impulse_responses(dsgevar, jlddata["data"], :cholesky, 1; use_intercept = true,
                              horizon = impulse_response_horizons(dsgevar))
    out_lr_h = impulse_responses(dsgevar, jlddata["data"], :choleskyLR, 1; use_intercept = true,
                              horizon = impulse_response_horizons(dsgevar))
    out_maxbc_h = impulse_responses(dsgevar, jlddata["data"], :maxBC, 1; use_intercept = true,
                              horizon = impulse_response_horizons(dsgevar))

    @test @test_matrix_approx_eq jlddata["exp_modal_cholesky_irf"] out
    @test @test_matrix_approx_eq jlddata["exp_modal_choleskyLR_irf"] out_lr
    @test @test_matrix_approx_eq jlddata["exp_modal_maxBC_irf"] out_maxbc
    @test @test_matrix_approx_eq jlddata["exp_modal_choleskyLR_irf"] out_lr2
    @test @test_matrix_approx_eq jlddata["exp_modal_maxBC_irf"] out_maxbc2
    @test @test_matrix_approx_eq jlddata["exp_modal_cholesky_irf"] -out_flip
    @test @test_matrix_approx_eq jlddata["exp_modal_choleskyLR_irf"] -out_lr_flip
    @test @test_matrix_approx_eq jlddata["exp_modal_maxBC_irf"] -out_maxbc_flip
    @test @test_matrix_approx_eq jlddata["exp_modal_cholesky_irf"] out_h
    @test @test_matrix_approx_eq jlddata["exp_modal_choleskyLR_irf"] out_lr_h
    @test @test_matrix_approx_eq jlddata["exp_modal_maxBC_irf"] out_maxbc_h
end

@testset "Impulse responses of VAR by using a DSGE as prior and to identify the rotation matrix" begin
    jlddata = load(joinpath(fp, "../../../reference/test_dsgevar_lambda_irfs.jld2"))
    m = Model1002("ss10", custom_settings =
                  Dict{Symbol,Setting}(:add_laborshare_measurement =>
                                       Setting(:add_laborshare_measurement, true)))
    m <= Setting(:impulse_response_horizons, 10)
    dsgevar = DSGEVAR(m, collect(keys(m.exogenous_shocks)), "ss11")
    DSGE.update!(dsgevar, λ = 1.)
    DSGE.update!(dsgevar, jlddata["modal_param"])

    Random.seed!(1793)
    out = impulse_responses(dsgevar, jlddata["data"])
    Random.seed!(1793)
    out_flip = impulse_responses(dsgevar, jlddata["data"]; flip_shocks = true)
    nobs = size(jlddata["data"], 1)
    lags = DSGE.get_lags(dsgevar)
    k = nobs * lags + 1
    XX = DSGE.lag_data(jlddata["data"], lags; use_intercept = true)
    X̂ = vcat(1, jlddata["data"][:, end], XX[end, 1+1:k - nobs])
    Random.seed!(1793)
    out_X̂ = impulse_responses(dsgevar, jlddata["data"], X̂)
    Random.seed!(1793)
    out_MM1 = impulse_responses(dsgevar, jlddata["data"]; MM = zeros(DSGE.n_observables(dsgevar),
                                                                     DSGE.n_shocks(dsgevar)))
    Random.seed!(1793)
    out_MM2 = impulse_responses(dsgevar, jlddata["data"]; MM = rand(DSGE.n_observables(dsgevar),
                                                                    DSGE.n_shocks(dsgevar)))
    Random.seed!(1793)
    out_draw = impulse_responses(dsgevar, jlddata["data"]; draw_shocks = true)


    @test @test_matrix_approx_eq jlddata["rotation_irf_by_shock"] out
    @test @test_matrix_approx_eq jlddata["flip_rotation_irf_by_shock"] out_flip
    @test @test_matrix_approx_eq out out_MM1
    @test !(out ≈ out_MM2)
    @test @test_matrix_approx_eq out out_X̂
    @test @test_matrix_approx_eq jlddata["rotation_irf_draw_shock"] out_draw
end

nothing
