writing_output = false
if VERSION < v"1.5"
    ver = "111"
else
    ver = "150"
end

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

if VERSION < v"1.5"
    @testset "Impulse responses to structural shocks identified by a DSGE-VAR" begin
        matdata1 = load(joinpath(fp, "../../../reference/test_irfdsge.jld2"))
        matdata2 = load(joinpath(fp, "../../../reference/test_rotations.jld2"))
        ŷ = DSGE.impulse_responses(matdata1["TTT"], matdata1["RRR"], matdata1["zz"],
                                   zeros(size(matdata1["zz"], 1)), matdata1["mmm"],
                                   matdata1["impact"].^2, Int(matdata2["k"]),
                                   convert(Matrix{Float64}, matdata2["cct_sim"]'),
                                   matdata2["sig_sim"], Int(matdata2["qahead"]),
                                   vec(matdata2["XXpred"]); normalize_rotation = false,
                                   test_shocks =
                                   convert(Matrix{Float64}, matdata2["Shocks"]'))

        @test ŷ ≈ matdata2["yypred"]'

        Random.seed!(1793) # drawing shocks, so this shouldn't be the same
        ŷ1 = DSGE.impulse_responses(matdata1["TTT"], matdata1["RRR"], matdata1["zz"],
                                    zeros(size(matdata1["zz"], 1)), matdata1["mmm"],
                                    matdata1["impact"].^2, Int(matdata2["k"]),
                                    convert(Matrix{Float64}, matdata2["cct_sim"]'),
                                    matdata2["sig_sim"], Int(matdata2["qahead"]),
                                    vec(matdata2["XXpred"]), draw_shocks = true,
                                    normalize_rotation = false)
        Random.seed!(1793) # re-seeding should work
        ŷ2 = DSGE.impulse_responses(matdata1["TTT"], matdata1["RRR"], matdata1["zz"],
                                    zeros(size(matdata1["zz"], 1)), matdata1["mmm"],
                                    matdata1["impact"].^2, Int(matdata2["k"]),
                                    convert(Matrix{Float64}, matdata2["cct_sim"]'),
                                    matdata2["sig_sim"], Int(matdata2["qahead"]),
                                    vec(matdata2["XXpred"]), draw_shocks = true,
                                    normalize_rotation = false)
        Random.seed!(1793) # flipping shocks shouldn't yield anything different
        @info "The following warning is expected."
        ŷ3 = DSGE.impulse_responses(matdata1["TTT"], matdata1["RRR"], matdata1["zz"],
                                    zeros(size(matdata1["zz"], 1)), matdata1["mmm"],
                                    matdata1["impact"].^2, Int(matdata2["k"]),
                                    convert(Matrix{Float64}, matdata2["cct_sim"]'),
                                    matdata2["sig_sim"], Int(matdata2["qahead"]),
                                    vec(matdata2["XXpred"]), draw_shocks = true,
                                    flip_shocks = true, normalize_rotation = false)
        @test !(ŷ1 ≈ matdata2["yypred"]')
        @test ŷ1 ≈ ŷ2
        @test ŷ1 ≈ ŷ3
    end

    @testset "Impulse responses of a VAR using a DSGE as a prior" begin
        jlddata = load(joinpath(fp, "../../../reference/test_dsgevar_lambda_irfs.jld2"))
        m = Model1002("ss10", custom_settings =
                      Dict{Symbol,Setting}(:add_laborshare_measurement =>
                                           Setting(:add_laborshare_measurement, true),
                                           :add_NominalWageGrowth =>
                                           Setting(:add_NominalWageGrowth, true)))
        m <= Setting(:impulse_response_horizons, 10)
        dsgevar = DSGEVAR(m, collect(keys(m.exogenous_shocks)), "ss11")
        DSGE.update!(dsgevar, λ = 1.)
        DSGE.update!(dsgevar, jlddata["modal_param"])

        Random.seed!(1793)
        out = impulse_responses(dsgevar, jlddata["data"], :cholesky, 1)
        out_lr = impulse_responses(dsgevar, jlddata["data"], :choleskyLR, 1)
        out_maxbc = impulse_responses(dsgevar, jlddata["data"], :maxBC, 1)

        Random.seed!(1793)
        _ = impulse_responses(dsgevar, jlddata["data"], :cholesky, 1)
        out_lr2 = impulse_responses(dsgevar, jlddata["data"], :cholesky_long_run, 1)
        out_maxbc2 = impulse_responses(dsgevar, jlddata["data"], :maximum_business_cycle_variance,
                                       1)

        Random.seed!(1793)
        out_flip = impulse_responses(dsgevar, jlddata["data"], :cholesky, 1,
                                     flip_shocks = true,)
        out_lr_flip = impulse_responses(dsgevar, jlddata["data"], :choleskyLR, 1,
                                        flip_shocks = true)
        out_maxbc_flip = impulse_responses(dsgevar, jlddata["data"], :maxBC, 1,
                                           flip_shocks = true)

        Random.seed!(1793)
        out_h = impulse_responses(dsgevar, jlddata["data"], :cholesky, 1,
                                  horizon = impulse_response_horizons(dsgevar))
        out_lr_h = impulse_responses(dsgevar, jlddata["data"], :choleskyLR, 1,
                                     horizon = impulse_response_horizons(dsgevar))
        out_maxbc_h = impulse_responses(dsgevar, jlddata["data"], :maxBC, 1,
                                        horizon = impulse_response_horizons(dsgevar))

        if writing_output
            jldopen(joinpath(fp, "../../../reference/test_dsgevar_lambda_irfs_output_version=" * ver * ".jld2"),
                    true, true, true, IOStream) do file
                        write(file, "exp_modal_cholesky_irf", out)
                        write(file, "exp_modal_choleskyLR_irf", out_lr)
                        write(file, "exp_modal_maxBC_irf", out_maxbc)
                    end
        end
        jlddata = load(joinpath(fp, "../../../reference/test_dsgevar_lambda_irfs_output_version=" * ver * ".jld2"))

        @test @test_matrix_approx_eq jlddata["exp_modal_cholesky_irf"] out
        @test @test_matrix_approx_eq jlddata["exp_modal_choleskyLR_irf"] out_lr

        @test @test_matrix_approx_eq jlddata["exp_modal_choleskyLR_irf"] out_lr2

        @test @test_matrix_approx_eq jlddata["exp_modal_cholesky_irf"] -out_flip
        @test @test_matrix_approx_eq jlddata["exp_modal_choleskyLR_irf"] -out_lr_flip

        @test @test_matrix_approx_eq jlddata["exp_modal_cholesky_irf"] out_h
        @test @test_matrix_approx_eq jlddata["exp_modal_choleskyLR_irf"] out_lr_h

        # Test maxBC separately b/c these have a slightly different error bound that leads to errors on Julia 1.0 but not Julia 1.1

        if VERSION >= v"1.1"
            @test @test_matrix_approx_eq jlddata["exp_modal_maxBC_irf"] out_maxbc
            @test @test_matrix_approx_eq jlddata["exp_modal_maxBC_irf"] out_maxbc2
            @test @test_matrix_approx_eq jlddata["exp_modal_maxBC_irf"] out_maxbc_h
            @test @test_matrix_approx_eq jlddata["exp_modal_maxBC_irf"] -out_maxbc_flip
        else
            @test maximum(abs.(jlddata["exp_modal_maxBC_irf"] - out_maxbc)) < 6e-6
            @test maximum(abs.(jlddata["exp_modal_maxBC_irf"] - out_maxbc2)) < 6e-6
            @test maximum(abs.(jlddata["exp_modal_maxBC_irf"] - out_maxbc_h)) < 6e-6
            @test maximum(abs.(jlddata["exp_modal_maxBC_irf"] + out_maxbc_flip)) < 6e-6
        end
    end

@testset "Impulse responses of VAR by using a DSGE as prior and to identify the rotation matrix" begin
    jlddata = load(joinpath(fp, "../../../reference/test_dsgevar_lambda_irfs.jld2"))
    m = Model1002("ss10", custom_settings =
                  Dict{Symbol,Setting}(:add_laborshare_measurement =>
                                       Setting(:add_laborshare_measurement, true),
                                       :add_NominalWageGrowth =>
                                       Setting(:add_NominalWageGrowth, true)))
    m <= Setting(:impulse_response_horizons, 10)
    dsgevar = DSGEVAR(m, collect(keys(m.exogenous_shocks)), "ss11")
    DSGE.update!(dsgevar, λ = 1.)
    DSGE.update!(dsgevar, jlddata["modal_param"])

    Random.seed!(1793)
    out = impulse_responses(dsgevar, jlddata["data"], normalize_rotation = false)
    Random.seed!(1793)
    out_flip = impulse_responses(dsgevar, jlddata["data"]; flip_shocks = true, normalize_rotation = false)
    nobs = size(jlddata["data"], 1)
    lags = DSGE.get_lags(dsgevar)
    k = nobs * lags + 1
    XX = DSGE.lag_data(jlddata["data"], lags; use_intercept = true)
    X̂ = vcat(1, jlddata["data"][:, end], XX[end, 1+1:k - nobs])
    Random.seed!(1793)
    out_X̂ = impulse_responses(dsgevar, jlddata["data"], X̂, normalize_rotation = false)
    Random.seed!(1793)
    out_MM1 = impulse_responses(dsgevar, jlddata["data"]; MM = zeros(DSGE.n_observables(dsgevar),
                                                                     DSGE.n_shocks(dsgevar)), normalize_rotation = false)
    Random.seed!(1793)
    out_MM2 = impulse_responses(dsgevar, jlddata["data"]; MM = rand(DSGE.n_observables(dsgevar),
                                                                    DSGE.n_shocks(dsgevar)), normalize_rotation = false)
    Random.seed!(1793)
    out_draw = impulse_responses(dsgevar, jlddata["data"]; draw_shocks = true, normalize_rotation = false)

    Random.seed!(1793) # now use deviations
    out_dev = impulse_responses(dsgevar, jlddata["data"], deviations = true, normalize_rotation = false)

    Random.seed!(1793) # flip shock
    out_dev_flip = impulse_responses(dsgevar, jlddata["data"], deviations = true,
                                     flip_shocks = true, normalize_rotation = false)

    Random.seed!(1793) # draw shock
    out_dev_draw = impulse_responses(dsgevar, jlddata["data"], deviations = true,
                                     draw_shocks = true, normalize_rotation = false)

    if writing_output
        jldopen(joinpath(fp, "../../../reference/test_dsgevar_lambda_irfs_output_rotation_version=" * ver * ".jld2"),
                true, true, true, IOStream) do file
            write(file, "rotation_irf_by_shock", out)
            write(file, "flip_rotation_irf_by_shock", out_flip)
            write(file, "rotation_irf_draw_shock", out_draw)
            write(file, "deviations_rotation_irf_by_shock", out_dev)
            write(file, "deviations_rotation_irf_draw_shock", out_dev_draw)
        end
    end

    jlddata = load(joinpath(fp, "../../../reference/test_dsgevar_lambda_irfs_output_rotation_version=" * ver * ".jld2"))

    @test @test_matrix_approx_eq jlddata["rotation_irf_by_shock"] out
    @test @test_matrix_approx_eq jlddata["flip_rotation_irf_by_shock"] out_flip
    @test @test_matrix_approx_eq out out_MM1
    @test !(out ≈ out_MM2)
    @test @test_matrix_approx_eq out out_X̂
    @test @test_matrix_approx_eq jlddata["rotation_irf_draw_shock"] out_draw
    @test @test_matrix_approx_eq out_dev jlddata["deviations_rotation_irf_by_shock"]
    @test @test_matrix_approx_eq out_dev_draw jlddata["deviations_rotation_irf_draw_shock"]
    @test @test_matrix_approx_eq out_dev -out_dev_flip
end

@testset "Impulse responses of a VAR approximation to a DSGE (or λ = ∞)" begin
    m = AnSchorfheide()
    m <= Setting(:impulse_response_horizons, 10)
    Random.seed!(1793)
    observables = [:obs_gdp, :obs_nominalrate, :z_t]

    shocks = collect(keys(m.exogenous_shocks))
    fp = dirname(@__FILE__)
    dsgevar = DSGEVAR(m, shocks, "ss0")
    DSGE.update!(dsgevar, lags = 4, observables = observables, λ = 1.)


    out1 = impulse_responses(dsgevar, :cholesky, 1)
    out2 = impulse_responses(dsgevar, :choleskyLR, 1)
    out3 = impulse_responses(dsgevar, :maximum_business_cycle_variance, 1)
    out4 = impulse_responses(dsgevar, :cholesky_long_run, 1)
    out5 = impulse_responses(dsgevar, :maxBC, 1)
    out6 = impulse_responses(dsgevar, :cholesky, 1,
                             flip_shocks = true)
    out7 = impulse_responses(dsgevar, :choleskyLR, 1,
                             flip_shocks = true)
    out8 = impulse_responses(dsgevar, :maxBC, 1,
                             flip_shocks = true)

    out9 = impulse_responses(dsgevar, :cholesky, 1,
                             use_intercept = true)
    out10 = impulse_responses(dsgevar, :choleskyLR, 1,
                              use_intercept = true)
    out11 = impulse_responses(dsgevar, :maxBC, 1,
                              use_intercept = true)
    out12 = impulse_responses(dsgevar, :cholesky, 1,
                              use_intercept = true, flip_shocks = true)
    out13 = impulse_responses(dsgevar, :choleskyLR, 1,
                              use_intercept = true, flip_shocks = true)
    out14 = impulse_responses(dsgevar, :maxBC, 1,
                              use_intercept = true, flip_shocks = true)

    jlddata = load(joinpath(fp, "../../../reference/var_approx_dsge_irfs.jld2"))

    @test @test_matrix_approx_eq jlddata["exp_cholesky"][:, :, 1] out1
    @test @test_matrix_approx_eq jlddata["exp_choleskyLR"][:, :, 1] out2
    @test @test_matrix_approx_eq jlddata["exp_maxBC"][:, :, 1] out3
    @test @test_matrix_approx_eq jlddata["exp_choleskyLR"][:, :, 1] out4
    @test @test_matrix_approx_eq jlddata["exp_maxBC"][:, :, 1] out5
    @test @test_matrix_approx_eq jlddata["exp_cholesky"][:, :, 1] -out6
    @test @test_matrix_approx_eq jlddata["exp_choleskyLR"][:, :, 1] -out7
    @test @test_matrix_approx_eq jlddata["exp_maxBC"][:, :, 1] -out8

    @test @test_matrix_approx_eq jlddata["exp_cholesky_int"][:, :, 1] out9
    @test @test_matrix_approx_eq jlddata["exp_choleskyLR_int"][:, :, 1] out10
    @test @test_matrix_approx_eq jlddata["exp_maxBC_int"][:, :, 1] out11
    @test @test_matrix_approx_eq jlddata["exp_cholesky_int"][:, :, 1] -out12
    @test @test_matrix_approx_eq jlddata["exp_choleskyLR_int"][:, :, 1] -out13
    @test @test_matrix_approx_eq jlddata["exp_maxBC_int"][:, :, 1] -out14
end

end
nothing
