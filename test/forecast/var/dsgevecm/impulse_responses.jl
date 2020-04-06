fp = dirname(@__FILE__)

@testset "Impulse responses to structural shocks identified by a DSGE-VECM" begin
    matdata = load(joinpath(fp, "../../../reference/dsgevecm_lambda_irfs.jld2"))
    ŷ = DSGE.impulse_responses(matdata["TTT"], matdata["RRR"], matdata["ZZ"],
                               vec(matdata["DD"]), matdata["MM"],
                               matdata["QQ"], Int(matdata["k"]), Int(matdata["nvar"]),
                               Int(matdata["coint"]), matdata["cct_sim"], matdata["sig_sim"],
                               matdata["cointvec"], vec(matdata["XXpred"]),
                               Int(matdata["qahead"]);
                               test_shocks = matdata["Shocks"])

    @test ŷ ≈ matdata["yypred"]

    Random.seed!(1793) # drawing shocks, so this shouldn't be the same
    ŷ1 = DSGE.impulse_responses(matdata["TTT"], matdata["RRR"], matdata["ZZ"],
                                vec(matdata["DD"]), matdata["MM"],
                                matdata["QQ"], Int(matdata["k"]), Int(matdata["nvar"]),
                                Int(matdata["coint"]), matdata["cct_sim"], matdata["sig_sim"],
                                matdata["cointvec"], vec(matdata["XXpred"]),
                                Int(matdata["qahead"]); draw_shocks = true)
    Random.seed!(1793) # test re-seeding works
    ŷ2 = DSGE.impulse_responses(matdata["TTT"], matdata["RRR"], matdata["ZZ"],
                                vec(matdata["DD"]), matdata["MM"],
                                matdata["QQ"], Int(matdata["k"]), Int(matdata["nvar"]),
                                Int(matdata["coint"]), matdata["cct_sim"], matdata["sig_sim"],
                                matdata["cointvec"], vec(matdata["XXpred"]),
                                Int(matdata["qahead"]); draw_shocks = true)
    Random.seed!(1793) # flipping shocks should yield the same thing
@info "The following warning is expected."
    ŷ3 = DSGE.impulse_responses(matdata["TTT"], matdata["RRR"], matdata["ZZ"],
                                vec(matdata["DD"]), matdata["MM"],
                                matdata["QQ"], Int(matdata["k"]), Int(matdata["nvar"]),
                                Int(matdata["coint"]), matdata["cct_sim"], matdata["sig_sim"],
                                matdata["cointvec"], vec(matdata["XXpred"]),
                                Int(matdata["qahead"]); draw_shocks = true,
                                flip_shocks = true)

    @test !(ŷ1 ≈ matdata["yypred"])
    @test ŷ1 ≈ ŷ2
    @test !(ŷ1 ≈ ŷ3)
end

# @testset "Impulse responses of a VAR using a DSGE as a prior" begin
#     jlddata = load(joinpath(fp, "../../../reference/test_dsgevar_lambda_irfs.jld2"))
#     m = Model1002("ss10", custom_settings =
#                   Dict{Symbol,Setting}(:add_laborshare_measurement =>
#                                        Setting(:add_laborshare_measurement, true)))
#     m <= Setting(:impulse_response_horizons, 10)
#     dsgevar = DSGEVAR(m, collect(keys(m.exogenous_shocks)), "ss11")
#     DSGE.update!(dsgevar, λ = 1.)
#     DSGE.update!(dsgevar, jlddata["modal_param"])

#     Random.seed!(1793)
#     out = impulse_responses(dsgevar, jlddata["data"], :cholesky, 1)
#     out_lr = impulse_responses(dsgevar, jlddata["data"], :choleskyLR, 1)
#     out_maxbc = impulse_responses(dsgevar, jlddata["data"], :maxBC, 1)

#     Random.seed!(1793)
#     _ = impulse_responses(dsgevar, jlddata["data"], :cholesky, 1)
#     out_lr2 = impulse_responses(dsgevar, jlddata["data"], :cholesky_long_run, 1)
#     out_maxbc2 = impulse_responses(dsgevar, jlddata["data"], :maximum_business_cycle_variance,
#                                    1)

#     Random.seed!(1793)
#     out_flip = impulse_responses(dsgevar, jlddata["data"], :cholesky, 1,
#                                  flip_shocks = true,)
#     out_lr_flip = impulse_responses(dsgevar, jlddata["data"], :choleskyLR, 1,
#                                     flip_shocks = true)
#     out_maxbc_flip = impulse_responses(dsgevar, jlddata["data"], :maxBC, 1,
#                                        flip_shocks = true)

#     Random.seed!(1793)
#     out_h = impulse_responses(dsgevar, jlddata["data"], :cholesky, 1,
#                               horizon = impulse_response_horizons(dsgevar))
#     out_lr_h = impulse_responses(dsgevar, jlddata["data"], :choleskyLR, 1,
#                               horizon = impulse_response_horizons(dsgevar))
#     out_maxbc_h = impulse_responses(dsgevar, jlddata["data"], :maxBC, 1,
#                               horizon = impulse_response_horizons(dsgevar))

#     @test @test_matrix_approx_eq jlddata["exp_modal_cholesky_irf"] out
#     @test @test_matrix_approx_eq jlddata["exp_modal_choleskyLR_irf"] out_lr

#     @test @test_matrix_approx_eq jlddata["exp_modal_choleskyLR_irf"] out_lr2

#     @test @test_matrix_approx_eq jlddata["exp_modal_cholesky_irf"] -out_flip
#     @test @test_matrix_approx_eq jlddata["exp_modal_choleskyLR_irf"] -out_lr_flip

#     @test @test_matrix_approx_eq jlddata["exp_modal_cholesky_irf"] out_h
#     @test @test_matrix_approx_eq jlddata["exp_modal_choleskyLR_irf"] out_lr_h

#     # Test maxBC separately b/c these have a slightly different error bound that leads to errors on Julia 1.0 but not Julia 1.1

#     if VERSION >= v"1.1"
#         @test @test_matrix_approx_eq jlddata["exp_modal_maxBC_irf"] out_maxbc
#         @test @test_matrix_approx_eq jlddata["exp_modal_maxBC_irf"] out_maxbc2
#         @test @test_matrix_approx_eq jlddata["exp_modal_maxBC_irf"] out_maxbc_h
#         @test @test_matrix_approx_eq jlddata["exp_modal_maxBC_irf"] -out_maxbc_flip
#     else
#         @test maximum(abs.(jlddata["exp_modal_maxBC_irf"] - out_maxbc)) < 5e-6
#         @test maximum(abs.(jlddata["exp_modal_maxBC_irf"] - out_maxbc2)) < 5e-6
#         @test maximum(abs.(jlddata["exp_modal_maxBC_irf"] - out_maxbc_h)) < 5e-6
#         @test maximum(abs.(jlddata["exp_modal_maxBC_irf"] + out_maxbc_flip)) < 5e-6
#     end
# end

# @testset "Impulse responses of VAR by using a DSGE as prior and to identify the rotation matrix" begin
#     jlddata = load(joinpath(fp, "../../../reference/test_dsgevar_lambda_irfs.jld2"))
#     m = Model1002("ss10", custom_settings =
#                   Dict{Symbol,Setting}(:add_laborshare_measurement =>
#                                        Setting(:add_laborshare_measurement, true)))
#     m <= Setting(:impulse_response_horizons, 10)
#     dsgevar = DSGEVAR(m, collect(keys(m.exogenous_shocks)), "ss11")
#     DSGE.update!(dsgevar, λ = 1.)
#     DSGE.update!(dsgevar, jlddata["modal_param"])

#     Random.seed!(1793)
#     out = impulse_responses(dsgevar, jlddata["data"])
#     Random.seed!(1793)
#     out_flip = impulse_responses(dsgevar, jlddata["data"]; flip_shocks = true)
#     nobs = size(jlddata["data"], 1)
#     lags = DSGE.get_lags(dsgevar)
#     k = nobs * lags + 1
#     XX = DSGE.lag_data(jlddata["data"], lags; use_intercept = true)
#     X̂ = vcat(1, jlddata["data"][:, end], XX[end, 1+1:k - nobs])
#     Random.seed!(1793)
#     out_X̂ = impulse_responses(dsgevar, jlddata["data"], X̂)
#     Random.seed!(1793)
#     out_MM1 = impulse_responses(dsgevar, jlddata["data"]; MM = zeros(DSGE.n_observables(dsgevar),
#                                                                      DSGE.n_shocks(dsgevar)))
#     Random.seed!(1793)
#     out_MM2 = impulse_responses(dsgevar, jlddata["data"]; MM = rand(DSGE.n_observables(dsgevar),
#                                                                     DSGE.n_shocks(dsgevar)))
#     Random.seed!(1793)
#     out_draw = impulse_responses(dsgevar, jlddata["data"]; draw_shocks = true)


#     @test @test_matrix_approx_eq jlddata["rotation_irf_by_shock"] out
#     @test @test_matrix_approx_eq jlddata["flip_rotation_irf_by_shock"] out_flip
#     @test @test_matrix_approx_eq out out_MM1
#     @test !(out ≈ out_MM2)
#     @test @test_matrix_approx_eq out out_X̂
#     @test @test_matrix_approx_eq jlddata["rotation_irf_draw_shock"] out_draw
# end

# @testset "Impulse responses of a VAR approximation to a DSGE (or λ = ∞)" begin
#     m = AnSchorfheide()
#     m <= Setting(:impulse_response_horizons, 10)
#     Random.seed!(1793)
#     observables = [:obs_gdp, :obs_nominalrate, :z_t]

#     shocks = collect(keys(m.exogenous_shocks))
#     fp = dirname(@__FILE__)
#     jlddata = load(joinpath(fp, "../../../reference/var_approx_dsge_irfs.jld2"))
#     dsgevar = DSGEVAR(m, shocks, "ss0")
#     DSGE.update!(dsgevar, lags = 4, observables = observables, λ = 1.)


#     out1 = impulse_responses(dsgevar, :cholesky, 1)
#     out2 = impulse_responses(dsgevar, :choleskyLR, 1)
#     out3 = impulse_responses(dsgevar, :maximum_business_cycle_variance, 1)
#     out4 = impulse_responses(dsgevar, :cholesky_long_run, 1)
#     out5 = impulse_responses(dsgevar, :maxBC, 1)
#     out6 = impulse_responses(dsgevar, :cholesky, 1,
#                              flip_shocks = true)
#     out7 = impulse_responses(dsgevar, :choleskyLR, 1,
#                              flip_shocks = true)
#     out8 = impulse_responses(dsgevar, :maxBC, 1,
#                              flip_shocks = true)

#     out9 = impulse_responses(dsgevar, :cholesky, 1,
#                              use_intercept = true)
#     out10 = impulse_responses(dsgevar, :choleskyLR, 1,
#                               use_intercept = true)
#     out11 = impulse_responses(dsgevar, :maxBC, 1,
#                               use_intercept = true)
#     out12 = impulse_responses(dsgevar, :cholesky, 1,
#                               use_intercept = true, flip_shocks = true)
#     out13 = impulse_responses(dsgevar, :choleskyLR, 1,
#                               use_intercept = true, flip_shocks = true)
#     out14 = impulse_responses(dsgevar, :maxBC, 1,
#                               use_intercept = true, flip_shocks = true)


#     @test @test_matrix_approx_eq jlddata["exp_cholesky"][:, :, 1] out1
#     @test @test_matrix_approx_eq jlddata["exp_choleskyLR"][:, :, 1] out2
#     @test @test_matrix_approx_eq jlddata["exp_maxBC"][:, :, 1] out3
#     @test @test_matrix_approx_eq jlddata["exp_choleskyLR"][:, :, 1] out4
#     @test @test_matrix_approx_eq jlddata["exp_maxBC"][:, :, 1] out5
#     @test @test_matrix_approx_eq jlddata["exp_cholesky"][:, :, 1] -out6
#     @test @test_matrix_approx_eq jlddata["exp_choleskyLR"][:, :, 1] -out7
#     @test @test_matrix_approx_eq jlddata["exp_maxBC"][:, :, 1] -out8

#     @test @test_matrix_approx_eq jlddata["exp_cholesky_int"][:, :, 1] out9
#     @test @test_matrix_approx_eq jlddata["exp_choleskyLR_int"][:, :, 1] out10
#     @test @test_matrix_approx_eq jlddata["exp_maxBC_int"][:, :, 1] out11
#     @test @test_matrix_approx_eq jlddata["exp_cholesky_int"][:, :, 1] -out12
#     @test @test_matrix_approx_eq jlddata["exp_choleskyLR_int"][:, :, 1] -out13
#     @test @test_matrix_approx_eq jlddata["exp_maxBC_int"][:, :, 1] -out14
# end

nothing
