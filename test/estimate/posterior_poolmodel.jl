# Note that this test asumes TPF properly works
using BenchmarkTools

pm = PoolModel("ss1")
pm <= Setting(:data_vintage, "190822")
filepath = dirname(@__FILE__)
pm <= Setting(:dataroot, "$(filepath)/../reference/")
data = df_to_matrix(pm, load_data(pm))

# This commented code produces the saved output
tuning = Dict(:r_star => 2., :c_init => 0.3, :target_accept_rate => 0.4,
              :resampling_method => :systematic, :n_mh_steps => 1,
              :n_particles => 1000, :n_presample_periods => 0,
              :allout => true)
pm <= Setting(:tuning, tuning, "tuning parameters for TPF")
pm[:ρ].value = 0.8
pm[:μ].value = 0.
pm[:σ].value = 1.
Random.seed!(1793)
Φpost, Ψpost, F_ϵpost, F_upost, F_λpost = compute_system(pm)
s_init = quantile.(Normal(), rand(F_λpost, tuning[:n_particles]))
# s_init = reshape(rand(F_λpost, tuning[:n_particles]), 1, 1000)
# s_init = [s_init; 1 .- s_init]
tpf_out, ~, ~ = tempered_particle_filter(data, Φpost, Ψpost, F_ϵpost, F_upost,
                                   s_init; tuning..., verbose = :none,
                                   fixed_sched = [1.], parallel = false, poolmodel = true)

@testset "Check likelihood and posterior calculations for dynamic weight" begin
    Random.seed!(1793)
    lh = likelihood(pm, data)
    @test lh == tpf_out

    Random.seed!(1793)
    x = map(α->α.value, pm.parameters)
    post_at_start = DSGE.posterior!(pm, x, data)
    @test post_at_start == tpf_out

    Random.seed!(1793)
    global y = x .+ [-.7; -300.; 5.]
    post_not_at_start = DSGE.posterior!(pm, y, data)
    ϵ = 0.0004
    @test abs(post_at_start - post_not_at_start) > ϵ
end

pm = PoolModel("ss1", weight_type = :static)
Random.seed!(1793)
static_out = sum(DSGE.filter_likelihood(pm, data))

@testset "Check likelihood and posterior calculations for static weight" begin
    Random.seed!(1793)
    lh = likelihood(pm, data)
    @test lh == static_out

    Random.seed!(1793)
    x = map(α->α.value, pm.parameters)
    post_at_start = DSGE.posterior!(pm, x, data)
    @test post_at_start == static_out

    Random.seed!(1793)
    global z = x .+ (1 .- x) ./ 2
    post_not_at_start = DSGE.posterior!(pm, z, data)
    ϵ = 1.
    @test abs(post_at_start - post_not_at_start) > ϵ
end

################
# Benchmarking #
################
# Flip to true to run; off by default.
run_benchmarks = false

if run_benchmarks
    bench_path = dirname(@__FILE__)
    bench_tuning = Dict(:r_star => 2., :c_init => 0.3, :target_accept_rate => 0.4,
                        :resampling_method => :systematic, :n_mh_steps => 1,
                        :n_particles => 1000, :n_presample_periods => 0,
                        :allout => true)

    # Dynamic pool: likelihood/posterior! 
    pm_dyn = PoolModel("ss1")
    pm_dyn <= Setting(:data_vintage, "190822")
    pm_dyn <= Setting(:dataroot, "$(bench_path)/../reference/")
    pm_dyn <= Setting(:tuning, bench_tuning, "tuning parameters for TPF")
    data_bench = df_to_matrix(pm_dyn, load_data(pm_dyn))
    x_dyn = map(α -> α.value, pm_dyn.parameters)

    # Static pool: closed-form score (cheap)
    pm_static = PoolModel("ss1", weight_type = :static)
    x_static = map(α -> α.value, pm_static.parameters)

    b_prior      = @benchmark prior($pm_dyn)
    b_lh_dyn     = @benchmark likelihood($pm_dyn, $data_bench)
    b_postx_dyn  = @benchmark DSGE.posterior!($pm_dyn, $x_dyn, $data_bench)
    b_lh_static  = @benchmark likelihood($pm_static, $data_bench)
    b_postx_stat = @benchmark DSGE.posterior!($pm_static, $x_static, $data_bench)

    println("\n===== estimate/posterior_poolmodel benchmark results =====")
    for (name, b) in [("prior              ", b_prior),
                      ("likelihood dynamic ", b_lh_dyn),
                      ("posterior! dynamic ", b_postx_dyn),
                      ("likelihood static  ", b_lh_static),
                      ("posterior! static  ", b_postx_stat)]
        println(name, "  time:   ", BenchmarkTools.prettytime(median(b).time),
                "   memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end

nothing
