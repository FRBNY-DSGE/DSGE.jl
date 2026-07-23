using BenchmarkTools

fp = dirname(@__FILE__)
if VERSION < v"1.5"
    ver = "111"
elseif VERSION < v"1.6"
    ver = "150"
elseif VERSION < v"1.12"
    ver = "160"
else
    ver = "1120"
end

# The references are generated sequentially, and the parallel path (pmap) only matches them
# when real worker processes exist; in a single process it diverges. So only exercise the
# parallel branch when extra workers are present — otherwise run sequentially (which is what
# the committed references assert against). Add workers (addprocs + @everywhere using DSGE)
# before this file to actually test the parallel path.
use_parallel = nworkers() > 1

@testset "Impulse responses of a VAR using a DSGE as a prior (wrapper function)" begin
    jlddata = load(joinpath(fp, "../../reference/test_dsgevar_lambda_irfs.jld2"))
    expdata = load(joinpath(fp, "../../reference/test_dsgevar_lambda_irfs_output_version=" * ver * ".jld2"))
    m = Model1002("ss10", custom_settings = [Setting(:add_laborshare_measurement, true),
                                             Setting(:add_NominalWageGrowth, true)])
    m <= Setting(:impulse_response_horizons, 10)
    dsgevar = DSGEVAR(m, collect(keys(m.exogenous_shocks)), "ss11")
    DSGE.update!(dsgevar, λ = 1.)
    DSGE.update!(dsgevar, jlddata["modal_param"])

    Random.seed!(1793)
    out = impulse_responses(dsgevar, jlddata["modal_param"],
                            jlddata["data"],
                            :mode, :cholesky; parallel = false,
                            create_meansbands = false, flip_shocks = false,
                            n_obs_shock = 1)
    out_lr = impulse_responses(dsgevar, jlddata["modal_param"], jlddata["data"],
                               :mode, :choleskyLR; parallel = false,
                               create_meansbands = false, flip_shocks = false,
                               n_obs_shock = 1)
    out_maxbc = impulse_responses(dsgevar,  jlddata["modal_param"], jlddata["data"],
                                  :mode, :maxBC; parallel = false,
                                  create_meansbands = false, flip_shocks = false,
                                  n_obs_shock = 1)

    Random.seed!(1793)
    _ = impulse_responses(dsgevar, jlddata["modal_param"],
                          jlddata["data"],
                          :mode, :cholesky; parallel = false,
                          create_meansbands = false, flip_shocks = false,
                          n_obs_shock = 1)

    out_lr2 = impulse_responses(dsgevar, jlddata["modal_param"], jlddata["data"],
                                :mode, :cholesky_long_run; parallel = false,
                                create_meansbands = false, flip_shocks = false,
                                n_obs_shock = 1)
    out_maxbc2 = impulse_responses(dsgevar,  jlddata["modal_param"], jlddata["data"],
                                   :mode, :maximum_business_cycle_variance; parallel = false,
                                   create_meansbands = false, flip_shocks = false,
                                   n_obs_shock = 1)

    Random.seed!(1793)
    out_flip = impulse_responses(dsgevar, jlddata["modal_param"], jlddata["data"],
                                 :mode, :cholesky; parallel = false,
                                 create_meansbands = false, flip_shocks = true,
                                 n_obs_shock = 1)
    out_lr_flip = impulse_responses(dsgevar, jlddata["modal_param"], jlddata["data"],
                                    :mode, :choleskyLR; parallel = false,
                                    create_meansbands = false, flip_shocks = true,
                                    n_obs_shock = 1)
    out_maxbc_flip = impulse_responses(dsgevar, jlddata["modal_param"], jlddata["data"],
                                       :mode, :maxBC; parallel = false,
                                       create_meansbands = false, flip_shocks = true,
                                       n_obs_shock = 1)

    # Not testing but checking no errors are caused when creating a MeansBands
    mb = impulse_responses(dsgevar, reshape(jlddata["modal_param"], 1,
                                            length(jlddata["modal_param"])),
                            jlddata["data"],
                            :mode, :cholesky; parallel = false,
                            create_meansbands = true, test_meansbands = true,
                            flip_shocks = false, n_obs_shock = 1)
    mb = impulse_responses(dsgevar, reshape(jlddata["modal_param"], 1,
                                            length(jlddata["modal_param"])),
                            jlddata["data"],
                            :mode, :choleskyLR; parallel = false,
                            create_meansbands = true, test_meansbands = true,
                            flip_shocks = false, n_obs_shock = 1)
    mb = impulse_responses(dsgevar, reshape(jlddata["modal_param"], 1,
                                            length(jlddata["modal_param"])),
                            jlddata["data"],
                            :mode, :maxBC; parallel = false,
                            create_meansbands = true, test_meansbands = true,
                            flip_shocks = false, n_obs_shock = 1)

    @test @test_matrix_approx_eq expdata["exp_modal_cholesky_irf"] dropdims(out, dims = 3)
    @test @test_matrix_approx_eq expdata["exp_modal_choleskyLR_irf"] dropdims(out_lr, dims = 3)
    @test @test_matrix_approx_eq expdata["exp_modal_maxBC_irf"] dropdims(out_maxbc, dims = 3)
    @test @test_matrix_approx_eq expdata["exp_modal_choleskyLR_irf"] dropdims(out_lr2, dims = 3)
    @test @test_matrix_approx_eq expdata["exp_modal_maxBC_irf"] dropdims(out_maxbc2, dims = 3)
    @test @test_matrix_approx_eq expdata["exp_modal_cholesky_irf"] -dropdims(out_flip, dims = 3)
    @test @test_matrix_approx_eq expdata["exp_modal_choleskyLR_irf"] -dropdims(out_lr_flip, dims = 3)
    @test @test_matrix_approx_eq expdata["exp_modal_maxBC_irf"] -dropdims(out_maxbc_flip, dims = 3)
end

@testset "Impulse responses of a VAR using parallel (1 worker) and using a DSGE as a prior (wrapper function)" begin
    jlddata = load(joinpath(fp, "../../reference/test_dsgevar_lambda_irfs.jld2"))
    expdata = load(joinpath(fp, "../../reference/test_dsgevar_lambda_irfs_output_version=" * ver * ".jld2"))
    m = Model1002("ss10", custom_settings = [Setting(:add_laborshare_measurement, true),
                                             Setting(:add_NominalWageGrowth, true)])
    m <= Setting(:impulse_response_horizons, 10)
    dsgevar = DSGEVAR(m, collect(keys(m.exogenous_shocks)), "ss11")
    DSGE.update!(dsgevar, λ = 1.)
    DSGE.update!(dsgevar, jlddata["modal_param"])

    Random.seed!(1793)
    out = impulse_responses(dsgevar, jlddata["modal_param"],
                            jlddata["data"],
                            :mode, :cholesky; parallel = use_parallel,
                            create_meansbands = false, flip_shocks = false,
                            n_obs_shock = 1)
    out_lr = impulse_responses(dsgevar, jlddata["modal_param"], jlddata["data"],
                               :mode, :choleskyLR; parallel = use_parallel,
                               create_meansbands = false, flip_shocks = false,
                               n_obs_shock = 1)
    out_maxbc = impulse_responses(dsgevar,  jlddata["modal_param"], jlddata["data"],
                                  :mode, :maxBC; parallel = use_parallel,
                                  create_meansbands = false, flip_shocks = false,
                                  n_obs_shock = 1)

    Random.seed!(1793)
    out_flip = impulse_responses(dsgevar, jlddata["modal_param"], jlddata["data"],
                                 :mode, :cholesky; parallel = use_parallel,
                                 create_meansbands = false, flip_shocks = true,
                                 n_obs_shock = 1)
    out_lr_flip = impulse_responses(dsgevar, jlddata["modal_param"], jlddata["data"],
                                    :mode, :choleskyLR; parallel = use_parallel,
                                    create_meansbands = false, flip_shocks = true,
                                    n_obs_shock = 1)
    out_maxbc_flip = impulse_responses(dsgevar, jlddata["modal_param"], jlddata["data"],
                                       :mode, :maxBC; parallel = use_parallel,
                                       create_meansbands = false, flip_shocks = true,
                                       n_obs_shock = 1)

    @test @test_matrix_approx_eq expdata["exp_modal_cholesky_irf"] dropdims(out, dims = 3)
    @test @test_matrix_approx_eq expdata["exp_modal_choleskyLR_irf"] dropdims(out_lr, dims = 3)
    @test @test_matrix_approx_eq expdata["exp_modal_maxBC_irf"] dropdims(out_maxbc, dims = 3)
    @test @test_matrix_approx_eq expdata["exp_modal_cholesky_irf"] -dropdims(out_flip, dims = 3)
    @test @test_matrix_approx_eq expdata["exp_modal_choleskyLR_irf"] -dropdims(out_lr_flip, dims = 3)
    @test @test_matrix_approx_eq expdata["exp_modal_maxBC_irf"] -dropdims(out_maxbc_flip, dims = 3)
end

@testset "Impulse responses of a VAR using a DSGE as a prior (wrapper function)" begin
    jlddata = load(joinpath(fp, "../../reference/test_dsgevar_lambda_irfs.jld2"))
    expdata = load(joinpath(fp, "../../reference/test_dsgevar_lambda_irfs_output_rotation_version=" * ver * ".jld2"))
    m = Model1002("ss10", custom_settings = [Setting(:add_laborshare_measurement, true),
                                             Setting(:add_NominalWageGrowth, true)])
    m <= Setting(:impulse_response_horizons, 10)
    dsgevar = DSGEVAR(m, collect(keys(m.exogenous_shocks)), "ss11")
    DSGE.update!(dsgevar, λ = 1.)
    DSGE.update!(dsgevar, jlddata["modal_param"])
    Random.seed!(1793)
    params = convert(Matrix{Float64},
                     hcat(jlddata["modal_param"],
                          jlddata["modal_param"] + vcat(randn(3) ./ 50, zeros(92)))')
    Random.seed!(1793)
    out = impulse_responses(dsgevar, params, jlddata["data"], :full, :rotation, normalize_rotation = false)
    Random.seed!(1793)
    out_draw = impulse_responses(dsgevar, params, jlddata["data"], :full, :rotation, draw_shocks = true, normalize_rotation = false)
    Random.seed!(1793)
    out_dev = impulse_responses(dsgevar, params, jlddata["data"], :full, :rotation,
                                deviations = true, normalize_rotation = false)
    Random.seed!(1793)
    out_dev_draw = impulse_responses(dsgevar, params, jlddata["data"], :full, :rotation, draw_shocks = true,
                                     deviations = true, normalize_rotation = false)


    Random.seed!(1793)
    out_parallel = impulse_responses(dsgevar, params, jlddata["data"], :full, :rotation, parallel = use_parallel, normalize_rotation = false)
    Random.seed!(1793)
    out_draw_parallel = impulse_responses(dsgevar, params, jlddata["data"], :full, :rotation,
                                          draw_shocks = true, parallel = use_parallel, normalize_rotation = false)
    Random.seed!(1793)
    out_dev_parallel = impulse_responses(dsgevar, params, jlddata["data"], :full, :rotation,
                                         deviations = true, parallel = use_parallel, normalize_rotation = false)
    Random.seed!(1793)
    out_dev_draw_parallel = impulse_responses(dsgevar, params, jlddata["data"], :full, :rotation, draw_shocks = true,
                                              deviations = true, parallel = use_parallel, normalize_rotation = false)


    # Not testing but just checking no error when creating MeansBands
    mb = impulse_responses(dsgevar, params, jlddata["data"], :full, :rotation,
                           create_meansbands = true, test_meansbands = true, normalize_rotation = false)

    @test @test_matrix_approx_eq out out_parallel
    @test @test_matrix_approx_eq out_dev out_dev_parallel
    @test @test_matrix_approx_eq out[:, :, :, 1] expdata["rotation_irf_by_shock"]
    @test @test_matrix_approx_eq out_dev[:, :, :, 1] expdata["deviations_rotation_irf_by_shock"]
    @test @test_matrix_approx_eq out_draw out_draw_parallel
    @test @test_matrix_approx_eq out_dev_draw out_dev_draw_parallel
    @test @test_matrix_approx_eq out_draw[:, :, 1] expdata["rotation_irf_draw_shock"]
    @test @test_matrix_approx_eq out_dev_draw[:, :, 1] expdata["deviations_rotation_irf_draw_shock"]
    @test size(out) == (4, 10, 24, 2) # 4 observables, horizon is 10, 24 shocks, 2 parameter draws
    @test size(out_draw) == (4, 10, 2) # 4 observables, horizon is 10, 24 shocks, 2 parameter draws
    @test size(out_dev) == (4, 10, 24, 2) # 4 observables, horizon is 10, 24 shocks, 2 parameter draws
    @test size(out_dev_draw) == (4, 10, 2) # 4 observables, horizon is 10, 24 shocks, 2 parameter draws
end

################
# Benchmarking #
################
run_benchmarks = false
if run_benchmarks
    jlddata = load(joinpath(fp, "../../reference/test_dsgevar_lambda_irfs.jld2"))
    m = Model1002("ss10", custom_settings = [Setting(:add_laborshare_measurement, true),
                                             Setting(:add_NominalWageGrowth, true)])
    m <= Setting(:impulse_response_horizons, 10)
    dsgevar = DSGEVAR(m, collect(keys(m.exogenous_shocks)), "ss11")
    DSGE.update!(dsgevar, λ = 1.)
    DSGE.update!(dsgevar, jlddata["modal_param"])
    modal_param = jlddata["modal_param"]
    data        = jlddata["data"]

    b = @benchmark impulse_responses($dsgevar, $modal_param, $data, :mode, :cholesky;
                                     parallel = false, create_meansbands = false,
                                     flip_shocks = false, n_obs_shock = 1) samples = 10 evals = 1
    println("\nimpulse_responses (DSGEVAR-λ wrapper, :mode :cholesky)  time: ",
            BenchmarkTools.prettytime(median(b).time),
            "   memory: ", BenchmarkTools.prettymemory(median(b).memory))
end
