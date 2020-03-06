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

@testset "Impulse responses of a VAR using a DSGE as a prior (wrapper function)" begin
    jlddata = load(joinpath(fp, "../../../reference/test_dsgevar_lambda_irfs.jld2"))
    m = Model1002("ss10", custom_settings =
                  Dict{Symbol,Setting}(:add_laborshare_measurement =>
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
                            n_obs_var = 1)
    out_lr = impulse_responses(dsgevar, jlddata["modal_param"], jlddata["data"],
                               :mode, :choleskyLR; parallel = false,
                               create_meansbands = false, flip_shocks = false,
                               n_obs_var = 1)
    out_maxbc = impulse_responses(dsgevar,  jlddata["modal_param"], jlddata["data"],
                                  :mode, :maxBC; parallel = false,
                                  create_meansbands = false, flip_shocks = false,
                                  n_obs_var = 1)

    Random.seed!(1793)
    _ = impulse_responses(dsgevar, jlddata["modal_param"],
                          jlddata["data"],
                          :mode, :cholesky; parallel = false,
                          create_meansbands = false, flip_shocks = false,
                          n_obs_var = 1)

    out_lr2 = impulse_responses(dsgevar, jlddata["modal_param"], jlddata["data"],
                                :mode, :cholesky_long_run; parallel = false,
                                create_meansbands = false, flip_shocks = false,
                                n_obs_var = 1)
    out_maxbc2 = impulse_responses(dsgevar,  jlddata["modal_param"], jlddata["data"],
                                   :mode, :maximum_business_cycle_variance; parallel = false,
                                   create_meansbands = false, flip_shocks = false,
                                   n_obs_var = 1)

    Random.seed!(1793)
    out_flip = impulse_responses(dsgevar, jlddata["modal_param"], jlddata["data"],
                                 :mode, :cholesky; parallel = false,
                                 create_meansbands = false, flip_shocks = true,
                                 n_obs_var = 1)
    out_lr_flip = impulse_responses(dsgevar, jlddata["modal_param"], jlddata["data"],
                                    :mode, :choleskyLR; parallel = false,
                                    create_meansbands = false, flip_shocks = true,
                                    n_obs_var = 1)
    out_maxbc_flip = impulse_responses(dsgevar, jlddata["modal_param"], jlddata["data"],
                                       :mode, :maxBC; parallel = false,
                                       create_meansbands = false, flip_shocks = true,
                                       n_obs_var = 1)

    @test @test_matrix_approx_eq jlddata["exp_modal_cholesky_irf"] dropdims(out, dims = 3)
    @test @test_matrix_approx_eq jlddata["exp_modal_choleskyLR_irf"] dropdims(out_lr, dims = 3)
    @test @test_matrix_approx_eq jlddata["exp_modal_maxBC_irf"] dropdims(out_maxbc, dims = 3)
    @test @test_matrix_approx_eq jlddata["exp_modal_choleskyLR_irf"] dropdims(out_lr2, dims = 3)
    @test @test_matrix_approx_eq jlddata["exp_modal_maxBC_irf"] dropdims(out_maxbc2, dims = 3)
    @test @test_matrix_approx_eq jlddata["exp_modal_cholesky_irf"] -dropdims(out_flip, dims = 3)
    @test @test_matrix_approx_eq jlddata["exp_modal_choleskyLR_irf"] -dropdims(out_lr_flip, dims = 3)
    @test @test_matrix_approx_eq jlddata["exp_modal_maxBC_irf"] -dropdims(out_maxbc_flip, dims = 3)

end

@testset "Impulse responses of a VAR using parallel (1 worker) and using a DSGE as a prior (wrapper function)" begin
    jlddata = load(joinpath(fp, "../../../reference/test_dsgevar_lambda_irfs.jld2"))
    m = Model1002("ss10", custom_settings =
                  Dict{Symbol,Setting}(:add_laborshare_measurement =>
                                       Setting(:add_laborshare_measurement, true)))
    m <= Setting(:impulse_response_horizons, 10)
    dsgevar = DSGEVAR(m, collect(keys(m.exogenous_shocks)), "ss11")
    DSGE.update!(dsgevar, λ = 1.)
    DSGE.update!(dsgevar, jlddata["modal_param"])

    Random.seed!(1793)
    out = impulse_responses(dsgevar, jlddata["modal_param"],
                            jlddata["data"],
                            :mode, :cholesky; parallel = true,
                            create_meansbands = false, flip_shocks = false,
                            n_obs_var = 1)
    out_lr = impulse_responses(dsgevar, jlddata["modal_param"], jlddata["data"],
                               :mode, :choleskyLR; parallel = true,
                               create_meansbands = false, flip_shocks = false,
                               n_obs_var = 1)
    out_maxbc = impulse_responses(dsgevar,  jlddata["modal_param"], jlddata["data"],
                                  :mode, :maxBC; parallel = true,
                                  create_meansbands = false, flip_shocks = false,
                                  n_obs_var = 1)

    Random.seed!(1793)
    out_flip = impulse_responses(dsgevar, jlddata["modal_param"], jlddata["data"],
                                 :mode, :cholesky; parallel = true,
                                 create_meansbands = false, flip_shocks = true,
                                 n_obs_var = 1)
    out_lr_flip = impulse_responses(dsgevar, jlddata["modal_param"], jlddata["data"],
                                    :mode, :choleskyLR; parallel = true,
                                    create_meansbands = false, flip_shocks = true,
                                    n_obs_var = 1)
    out_maxbc_flip = impulse_responses(dsgevar, jlddata["modal_param"], jlddata["data"],
                                       :mode, :maxBC; parallel = true,
                                       create_meansbands = false, flip_shocks = true,
                                       n_obs_var = 1)

    @test @test_matrix_approx_eq jlddata["exp_modal_cholesky_irf"] dropdims(out, dims = 3)
    @test @test_matrix_approx_eq jlddata["exp_modal_choleskyLR_irf"] dropdims(out_lr, dims = 3)
    @test @test_matrix_approx_eq jlddata["exp_modal_maxBC_irf"] dropdims(out_maxbc, dims = 3)
    @test @test_matrix_approx_eq jlddata["exp_modal_cholesky_irf"] -dropdims(out_flip, dims = 3)
    @test @test_matrix_approx_eq jlddata["exp_modal_choleskyLR_irf"] -dropdims(out_lr_flip, dims = 3)
    @test @test_matrix_approx_eq jlddata["exp_modal_maxBC_irf"] -dropdims(out_maxbc_flip, dims = 3)
end


@testset "Impulse responses of a VAR using a DSGE as a prior (wrapper function)" begin
    jlddata = load(joinpath(fp, "../../../reference/test_dsgevar_lambda_irfs.jld2"))
    m = Model1002("ss10", custom_settings =
                  Dict{Symbol,Setting}(:add_laborshare_measurement =>
                                       Setting(:add_laborshare_measurement, true)))
    m <= Setting(:impulse_response_horizons, 10)
    dsgevar = DSGEVAR(m, collect(keys(m.exogenous_shocks)), "ss11")
    DSGE.update!(dsgevar, λ = 1.)
    DSGE.update!(dsgevar, jlddata["modal_param"])
    Random.seed!(1793)
    params = convert(Matrix{Float64},
                     hcat(jlddata["modal_param"],
                          jlddata["modal_param"] + vcat(randn(3) ./ 50, zeros(92)))')
    Random.seed!(1793)
    out = impulse_responses(dsgevar, params, jlddata["data"], :full, :rotation)
    Random.seed!(1793)
    out_draw = impulse_responses(dsgevar, params, jlddata["data"], :full, :rotation, draw_shocks = true)

    Random.seed!(1793)
    out_parallel = impulse_responses(dsgevar, params, jlddata["data"], :full, :rotation, parallel = true)
    Random.seed!(1793)
    out_draw_parallel = impulse_responses(dsgevar, params, jlddata["data"], :full, :rotation,
                                          draw_shocks = true, parallel = true)

    @test @test_matrix_approx_eq out out_parallel
    @test @test_matrix_approx_eq out_draw out_draw_parallel
    @test size(out) == (4, 10, 24, 2) # 4 observables, horizon is 10, 24 shocks, 2 parameter draws
    @test size(out_draw) == (4, 10, 2) # 4 observables, horizon is 10, 24 shocks, 2 parameter draws
end

nothing
