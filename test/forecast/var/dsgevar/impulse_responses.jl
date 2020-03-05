fp = dirname(@__FILE__)
@testset "DSGE impulse responses to a pre-specified impact matrix" begin
    matdata = load(joinpath(fp, "../../../reference/test_irfdsge.jld2"))
    irfout = impulse_responses(matdata["TTT"], matdata["RRR"], matdata["zz"],
                               zeros(size(matdata["zz"], 1)), matdata["mmm"],
                               matdata["impact"], 1; accumulate = true,
                               cum_inds = 1)
    nobs = size(matdata["zz"], 1)
    for i = 1:nobs
        @test vec(matdata["aairf"][:, nobs * (i - 1) + 1:nobs * i]) ≈
            vec(irfout[:, :, i])
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
                               Int(matdata2["qahead"]); accumulate = true,
                               cum_inds = 1, test_shocks =
                               convert(Matrix{Float64}, matdata2["Shocks"]'))

    @test ŷ ≈ matdata2["yypred"]
end

@testset "Impulse responses of a VAR using a DSGE as a prior" begin
    jlddata = load(joinpath(fp, "../../../reference/test_dsgevar_lambda_irfs.jld2"))
    m = Model1002("ss10", custom_settings = Dict{Symbol,Setting}(:add_laborshare_measurement =>
                                                                 Setting(:add_laborshare_measurement, true)))
    m <= Setting(:impulse_response_horizons, 10)
    dsgevar = DSGEVAR(m, collect(keys(m.exogenous_shocks)), "ss11")
    DSGE.update!(dsgevar, λ = 1.)
    DSGE.update!(dsgevar, jlddata["modal_param"])

    Random.seed!(1793)
    out = impulse_responses(dsgevar, jlddata["modal_param"],
                            jlddata["data"],
                            :mode, :cholesky; parallel = false,
                            create_meansbands = false, flip_shocks = false,
                            use_intercept = true, n_obs_var = 1)
    out_lr = impulse_responses(dsgevar, jlddata["modal_param"], jlddata["data"],
                               :mode, :choleskyLR; parallel = false,
                               create_meansbands = false, flip_shocks = false,
                               use_intercept = true, n_obs_var = 1)
    out_maxbc = impulse_responses(dsgevar,  jlddata["modal_param"], jlddata["data"],
                                  :mode, :maxBC; parallel = false,
                                  create_meansbands = false, flip_shocks = false,
                                  use_intercept = true, n_obs_var = 1)

    Random.seed!(1793)
    out_flip = impulse_responses(dsgevar, jlddata["modal_param"], jlddata["data"],
                                 :mode, :cholesky; parallel = false,
                                 create_meansbands = false, flip_shocks = true,
                                 use_intercept = true, n_obs_var = 1)
    out_lr_flip = impulse_responses(dsgevar, jlddata["modal_param"], jlddata["data"],
                                    :mode, :choleskyLR; parallel = false,
                                    create_meansbands = false, flip_shocks = true,
                                    use_intercept = true, n_obs_var = 1)
    out_maxbc_flip = impulse_responses(dsgevar, jlddata["modal_param"], jlddata["data"],
                                       :mode, :maxBC; parallel = false,
                                       create_meansbands = false, flip_shocks = true,
                                       use_intercept = true, n_obs_var = 1)

    @test @test_matrix_approx_eq jlddata["exp_modal_cholesky_irf"] dropdims(out, dims = 3)
    @test @test_matrix_approx_eq jlddata["exp_modal_choleskyLR_irf"] dropdims(out_lr, dims = 3)
    @test @test_matrix_approx_eq jlddata["exp_modal_maxBC_irf"] dropdims(out_maxbc, dims = 3)
    @test @test_matrix_approx_eq jlddata["exp_modal_cholesky_irf"] -dropdims(out_flip, dims = 3)
    @test @test_matrix_approx_eq jlddata["exp_modal_choleskyLR_irf"] -dropdims(out_lr_flip, dims = 3)
    @test @test_matrix_approx_eq jlddata["exp_modal_maxBC_irf"] -dropdims(out_maxbc_flip, dims = 3)
end

nothing
