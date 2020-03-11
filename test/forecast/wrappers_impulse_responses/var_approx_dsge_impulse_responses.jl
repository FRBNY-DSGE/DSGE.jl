# Test on DSGE
m = AnSchorfheide()
m <= Setting(:impulse_response_horizons, 10)
paras = map(x -> x.value, m.parameters)
Random.seed!(1793)
paras = convert(Matrix{Float64}, hcat(paras, paras + vcat(randn(2) ./ 100, zeros(14)))')
observables = [:obs_gdp, :obs_nominalrate, :z_t]

shocks = collect(keys(m.exogenous_shocks))
fp = dirname(@__FILE__)
jlddata = load(joinpath(fp, "../../reference/var_approx_dsge_irfs.jld2"))

@testset "Identified impulse responses from shocks to observables for a DSGE" begin
    out1 = impulse_responses(m, paras, :full, :cholesky, 4, observables, shocks, 1)
    out2 = impulse_responses(m, paras, :full, :choleskyLR, 4, observables, shocks, 1)
    out3 = impulse_responses(m, paras, :full, :maximum_business_cycle_variance, 4,
                             observables, shocks, 1)
    out4 = impulse_responses(m, paras, :full, :cholesky_long_run, 4, observables, shocks, 1)
    out5 = impulse_responses(m, paras, :full, :maxBC, 4, observables, shocks, 1)
    out6 = impulse_responses(m, paras, :full, :cholesky, 4, observables, shocks, 1,
                             flip_shocks = true)
    out7 = impulse_responses(m, paras, :full, :choleskyLR, 4, observables, shocks, 1,
                             flip_shocks = true)
    out8 = impulse_responses(m, paras, :full, :maxBC, 4, observables, shocks, 1,
                             flip_shocks = true)

    out9 = impulse_responses(m, paras, :full, :cholesky, 4, observables, shocks, 1,
                             use_intercept = true)
    out10 = impulse_responses(m, paras, :full, :choleskyLR, 4, observables, shocks, 1,
                              use_intercept = true)
    out11 = impulse_responses(m, paras, :full, :maxBC, 4, observables, shocks, 1,
                              use_intercept = true)
    out12 = impulse_responses(m, paras, :full, :cholesky, 4, observables, shocks, 1,
                              use_intercept = true, flip_shocks = true)
    out13 = impulse_responses(m, paras, :full, :choleskyLR, 4, observables, shocks, 1,
                              use_intercept = true, flip_shocks = true)
    out14 = impulse_responses(m, paras, :full, :maxBC, 4, observables, shocks, 1,
                              use_intercept = true, flip_shocks = true)

    out15 = impulse_responses(m, paras, :full, :cholesky, 4, observables, shocks, 1,
                              parallel = true)
    out16 = impulse_responses(m, paras, :full, :choleskyLR, 4, observables, shocks, 1,
                              parallel = true)
    out17 = impulse_responses(m, paras, :full, :maxBC, 4, observables, shocks, 1,
                              parallel = true)

    out18 = impulse_responses(m, paras[1, :], :mode, :cholesky, 4, observables, shocks, 1)
    out19 = impulse_responses(m, paras[1, :], :mode, :choleskyLR, 4, observables, shocks, 1)
    out20 = impulse_responses(m, paras[1, :], :mode, :maximum_business_cycle_variance, 4,
                              observables, shocks, 1)

    @test @test_matrix_approx_eq jlddata["exp_cholesky"] out1
    @test @test_matrix_approx_eq jlddata["exp_choleskyLR"] out2
    @test @test_matrix_approx_eq jlddata["exp_maxBC"] out3
    @test @test_matrix_approx_eq jlddata["exp_choleskyLR"] out4
    @test @test_matrix_approx_eq jlddata["exp_maxBC"] out5
    @test @test_matrix_approx_eq jlddata["exp_cholesky"] -out6
    @test @test_matrix_approx_eq jlddata["exp_choleskyLR"] -out7
    @test @test_matrix_approx_eq jlddata["exp_maxBC"] -out8

    @test @test_matrix_approx_eq jlddata["exp_cholesky_int"] out9
    @test @test_matrix_approx_eq jlddata["exp_choleskyLR_int"] out10
    @test @test_matrix_approx_eq jlddata["exp_maxBC_int"] out11
    @test @test_matrix_approx_eq jlddata["exp_cholesky_int"] -out12
    @test @test_matrix_approx_eq jlddata["exp_choleskyLR_int"] -out13
    @test @test_matrix_approx_eq jlddata["exp_maxBC_int"] -out14

    @test @test_matrix_approx_eq jlddata["exp_cholesky"] out15
    @test @test_matrix_approx_eq jlddata["exp_choleskyLR"] out16
    @test @test_matrix_approx_eq jlddata["exp_maxBC"] out17

    @test @test_matrix_approx_eq jlddata["exp_cholesky"][:, :, 1] out18
    @test @test_matrix_approx_eq jlddata["exp_choleskyLR"][:, :, 1] out19
    @test @test_matrix_approx_eq jlddata["exp_maxBC"][:, :, 1] out20
end

# Test on a DSGEVAR.
dsgevar = DSGEVAR(m, shocks, "ss0")
DSGE.update!(dsgevar, lags = 4, observables = observables, λ = 1.)


@testset "Impulse responses identified from shocks to observables for a DSGE when using a DSGEVAR" begin
    out1 = impulse_responses(dsgevar, paras, :full, :cholesky, 1)
    out2 = impulse_responses(dsgevar, paras, :full, :choleskyLR, 1)
    out3 = impulse_responses(dsgevar, paras, :full, :maximum_business_cycle_variance, 1)
    out4 = impulse_responses(dsgevar, paras, :full, :cholesky_long_run, 1)
    out5 = impulse_responses(dsgevar, paras, :full, :maxBC, 1)
    out6 = impulse_responses(dsgevar, paras, :full, :cholesky, 1,
                             flip_shocks = true)
    out7 = impulse_responses(dsgevar, paras, :full, :choleskyLR, 1,
                             flip_shocks = true)
    out8 = impulse_responses(dsgevar, paras, :full, :maxBC, 1,
                             flip_shocks = true)

    out9 = impulse_responses(dsgevar, paras, :full, :cholesky, 1,
                             use_intercept = true)
    out10 = impulse_responses(dsgevar, paras, :full, :choleskyLR, 1,
                              use_intercept = true)
    out11 = impulse_responses(dsgevar, paras, :full, :maxBC, 1,
                              use_intercept = true)
    out12 = impulse_responses(dsgevar, paras, :full, :cholesky, 1,
                              use_intercept = true, flip_shocks = true)
    out13 = impulse_responses(dsgevar, paras, :full, :choleskyLR, 1,
                              use_intercept = true, flip_shocks = true)
    out14 = impulse_responses(dsgevar, paras, :full, :maxBC, 1,
                              use_intercept = true, flip_shocks = true)

    out15 = impulse_responses(dsgevar, paras, :full, :cholesky, 1,
                              parallel = true)
    out16 = impulse_responses(dsgevar, paras, :full, :choleskyLR, 1,
                              parallel = true)
    out17 = impulse_responses(dsgevar, paras, :full, :maxBC, 1,
                              parallel = true)

    out18 = impulse_responses(dsgevar, paras[1, :], :mode, :cholesky, 1)
    out19 = impulse_responses(dsgevar, paras[1, :], :mode, :choleskyLR, 1)
    out20 = impulse_responses(dsgevar, paras[1, :], :mode, :maxBC, 1)

    out21 = impulse_responses(dsgevar, paras, :full, :cholesky, 4,
                              observables, shocks, 1)
    out22 = impulse_responses(dsgevar, paras, :full, :choleskyLR, 4,
                              observables, shocks, 1)
    out23 = impulse_responses(dsgevar, paras, :full, :maxBC, 4,
                              observables, shocks, 1)

    out24 = impulse_responses(dsgevar, paras[1, :], :mode, :cholesky, 4,
                              observables, shocks, 1)
    out25 = impulse_responses(dsgevar, paras[1, :], :mode, :choleskyLR, 4,
                              observables, shocks, 1)
    out26 = impulse_responses(dsgevar, paras[1, :], :mode, :maxBC, 4,
                              observables, shocks, 1)

    # Not testing but just checking for no errors when creating a MeansBands
    mb = impulse_responses(dsgevar, paras, :full, :cholesky, 1,
                           create_meansbands = true, test_meansbands = true)
    mb = impulse_responses(dsgevar, paras, :full, :choleskyLR, 1,
                           create_meansbands = true, test_meansbands = true)
    mb = impulse_responses(dsgevar, paras, :full, :maxBC, 1,
                           create_meansbands = true, test_meansbands = true)

    @test @test_matrix_approx_eq jlddata["exp_cholesky"] out1
    @test @test_matrix_approx_eq jlddata["exp_choleskyLR"] out2
    @test @test_matrix_approx_eq jlddata["exp_maxBC"] out3
    @test @test_matrix_approx_eq jlddata["exp_choleskyLR"] out4
    @test @test_matrix_approx_eq jlddata["exp_maxBC"] out5
    @test @test_matrix_approx_eq jlddata["exp_cholesky"] -out6
    @test @test_matrix_approx_eq jlddata["exp_choleskyLR"] -out7
    @test @test_matrix_approx_eq jlddata["exp_maxBC"] -out8

    @test @test_matrix_approx_eq jlddata["exp_cholesky_int"] out9
    @test @test_matrix_approx_eq jlddata["exp_choleskyLR_int"] out10
    @test @test_matrix_approx_eq jlddata["exp_maxBC_int"] out11
    @test @test_matrix_approx_eq jlddata["exp_cholesky_int"] -out12
    @test @test_matrix_approx_eq jlddata["exp_choleskyLR_int"] -out13
    @test @test_matrix_approx_eq jlddata["exp_maxBC_int"] -out14

    @test @test_matrix_approx_eq jlddata["exp_cholesky"] out15
    @test @test_matrix_approx_eq jlddata["exp_choleskyLR"] out16
    @test @test_matrix_approx_eq jlddata["exp_maxBC"] out17

    @test @test_matrix_approx_eq jlddata["exp_cholesky"][:, :, 1] out18
    @test @test_matrix_approx_eq jlddata["exp_choleskyLR"][:, :, 1] out19
    @test @test_matrix_approx_eq jlddata["exp_maxBC"][:, :, 1] out20

    @test @test_matrix_approx_eq jlddata["exp_cholesky"] out21
    @test @test_matrix_approx_eq jlddata["exp_choleskyLR"] out22
    @test @test_matrix_approx_eq jlddata["exp_maxBC"] out23

    @test @test_matrix_approx_eq jlddata["exp_cholesky"][:, :, 1] out24
    @test @test_matrix_approx_eq jlddata["exp_choleskyLR"][:, :, 1] out25
    @test @test_matrix_approx_eq jlddata["exp_maxBC"][:, :, 1] out26
end
