using BenchmarkTools
fp = dirname(@__FILE__)

@testset "VAR impulse responses" begin
    data = load(joinpath(fp, "../../reference/var_irf_out.jld2"))
    β = data["bhat_full"]
    Σ = data["sigmahat_full"]
    irfs1 = data["irfs_full_choleskyLR"]
    irfs2 = data["irfs_full_maxBC"]
    irfs3 = data["irfs_full_cholesky"]
    out1 = impulse_responses(β, Σ, 1, 20;
                             method = :choleskyLR,
                             flip_shocks = true)
    out2 = impulse_responses(β, Σ, 1, 20;
                             method = :maxBC,
                             flip_shocks = true)
    out3 = impulse_responses(β, Σ, 1, 20;
                             method = :cholesky,
                             flip_shocks = true)
    out4 = impulse_responses(β, Σ, 1, 20;
                             method = :cholesky_long_run,
                             flip_shocks = true)
    out5 = impulse_responses(β, Σ, 1, 20;
                             method = :maximum_business_cycle_variance,
                             flip_shocks = true)
    out6 = impulse_responses(β, Σ, 1, 20;
                             method = :choleskyLR,
                             flip_shocks = false)
    out7 = impulse_responses(β, Σ, 1, 20;
                             method = :maxBC,
                             flip_shocks = false)
    out8 = impulse_responses(β, Σ, 1, 20;
                             method = :cholesky,
                             flip_shocks = false)

    @test @test_matrix_approx_eq irfs1' out1
    @test @test_matrix_approx_eq irfs2' out2
    @test @test_matrix_approx_eq irfs3' out3
    @test @test_matrix_approx_eq irfs1' out4
    @test @test_matrix_approx_eq irfs2' out5
    @test @test_matrix_approx_eq irfs1[1, :] -out6[:, 1]
    @test @test_matrix_approx_eq irfs2[1, :] -out7[:, 1]
    @test @test_matrix_approx_eq irfs3[1, :] -out8[:, 1]
end


@testset "VECM impulse responses" begin
    data = load(joinpath(fp, "../../reference/dsgevecm_lambda_irfs.jld2"))
    expout = load(joinpath(fp, "../../reference/vecm_irfs.jld2"), "expout")
    β = data["cct_sim"]
    Σ = data["sig_sim"]
    coint_mat = data["cointvec"]
    out1 = impulse_responses(β, Σ, coint_mat, 1, 5;
                             method = :cholesky,
                             flip_shocks = true,
                             use_intercept = true)
    out2 = impulse_responses(β, Σ, coint_mat, 1, 5;
                             method = :cholesky,
                             flip_shocks = false,
                             use_intercept = true)
    out3 = impulse_responses(β[vcat(1:3, 5:32), :], Σ, coint_mat, 1, 5;
                             method = :cholesky, use_intercept = false,
                             flip_shocks = true)
    out4 = impulse_responses(β[vcat(1:3, 5:32), :], Σ, coint_mat, 1, 5;
                             method = :cholesky, use_intercept = false,
                             flip_shocks = false)


    @test_throws ErrorException impulse_responses(β, Σ, coint_mat, 1, 5;
                                                  method = :not_a_method,
                                                  use_intercept = true)

    @test @test_matrix_approx_eq expout out1
    @test @test_matrix_approx_eq expout -out2
    @test @test_matrix_approx_eq expout out3
    @test @test_matrix_approx_eq expout -out4
end

run_benchmarks = false
if run_benchmarks
    # Minimal top-level setup (testset-local vars are not visible here): reload the same
    # reference inputs and benchmark the primary VAR impulse_responses call (Cholesky).
    bench_data = load(joinpath(fp, "../../reference/var_irf_out.jld2"))
    bench_β = bench_data["bhat_full"]
    bench_Σ = bench_data["sigmahat_full"]

    b_var_irf = @benchmark impulse_responses($bench_β, $bench_Σ, 1, 20;
                                             method = :cholesky, flip_shocks = true)
    println("\n===== VAR impulse_responses benchmark results =====")
    println(rpad("impulse_responses (VAR, :cholesky)", 36), " time: ",
            rpad(BenchmarkTools.prettytime(median(b_var_irf).time), 12),
            "memory: ", BenchmarkTools.prettymemory(median(b_var_irf).memory))
end

nothing
