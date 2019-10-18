# Note that this test assumes TPF properly works
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
pm[:ρ].value = 0.5
pm[:μ].value = 0.
pm[:σ].value = 1.
Random.seed!(1793)
Φfilt, Ψfilt, F_ϵfilt, F_ufilt, F_λfilt = compute_system(pm)
s_init = reshape(rand(F_λfilt, tuning[:n_particles]), 1, 1000)
s_init = [s_init; 1 .- s_init]
tpf_out, ~, ~ = tempered_particle_filter(data, Φfilt, Ψfilt, F_ϵfilt, F_ufilt,
                                   s_init; tuning..., verbose = :none,
                                   fixed_sched = [1.], parallel = false, poolmodel = true)

Random.seed!(1793)
filt_tpf_out, ~, ~ = DSGE.filter(pm, data; tuning = get_setting(pm, :tuning))
Random.seed!(1793)
filt_lik_tpf_out = sum(DSGE.filter_likelihood(pm, data; tuning = get_setting(pm, :tuning)))

@testset "Check call to tempered particle filter without providing initial states" begin
    @test tpf_out == filt_tpf_out
    @test tpf_out == filt_lik_tpf_out
end

Random.seed!(1793)
s_init = reshape(rand(F_λfilt, tuning[:n_particles]), 1, 1000)
s_init = [s_init; 1 .- s_init] # this tpf output should be saved later
filt_tpf_out, ~, ~ = DSGE.filter(pm, data, s_init; tuning = get_setting(pm, :tuning))
Random.seed!(1793)
s_init = reshape(rand(F_λfilt, tuning[:n_particles]), 1, 1000)
s_init = [s_init; 1 .- s_init] # this tpf output should be saved later
filt_lik_tpf_out = sum(DSGE.filter_likelihood(pm, data, s_init; tuning = get_setting(pm, :tuning)))

@testset "Check call to tempered particle filter when providing initial states" begin
    @test tpf_out == filt_tpf_out
    @test tpf_out == filt_lik_tpf_out
end

pm = PoolModel("ss1", weight_type = :equal)
filepath = dirname(@__FILE__)
pm <= Setting(:dataroot, "$(filepath)/../reference/")
pm <= Setting(:data_vintage, "190822")
data = df_to_matrix(pm, load_data(pm))
filt_out, ~ = DSGE.filter(pm, data[:,1:2])
filtlik_out = sum(DSGE.filter_likelihood(pm, data[:,1:2]))
truefilt_out = sum(log.(sum(data[:,1:2] ./ 2, dims = 1)))
@testset "Check equal weights computes properly" begin
    @test filt_out == truefilt_out
    @test filtlik_out == truefilt_out
end

pm = PoolModel("ss1", weight_type = :static)
pm <= Setting(:data_vintage, "190822")
filepath = dirname(@__FILE__)
pm <= Setting(:dataroot, "$(filepath)/../reference/")
pm <= Setting(:data_vintage, "190822")
data = df_to_matrix(pm, load_data(pm))
s_init = [0.5, 0.5]
filt_out, ~  = DSGE.filter(pm, data[:,1:2], s_init)
filtlik_out = sum(DSGE.filter_likelihood(pm, data[:,1:2], s_init))
truefilt_out = sum(log.(sum(data[:,1:2] ./ 2, dims = 1)))
@testset "Check static weights computes properly" begin
    @test filt_out == truefilt_out
    @test filtlik_out == truefilt_out
end

@testset "Check warning and error conditions" begin
    pm = PoolModel("ss1")
    pm <= Setting(:data_vintage, "190822")
    filepath = dirname(@__FILE__)
    pm <= Setting(:dataroot, "$(filepath)/../reference/")
    tuning = get_setting(pm, :tuning)
    tuning[:n_particles] = 10
    pm <= Setting(:tuning, tuning)
    @test_throws ErrorException DSGE.filter(pm, data, vec([1. 1.; 1. 1.]))
    delete!(pm.settings, :tuning)
    @test_logs (:warn, "no tuning parameters provided; using default tempered particle filter values") DSGE.filter(pm, data)
    pm = PoolModel("ss1", weight_type = :bma)
    pm <= Setting(:data_vintage, "190822")
    filepath = dirname(@__FILE__)
    pm <= Setting(:dataroot, "$(filepath)/../reference/")
    @test_throws ErrorException DSGE.filter(pm, data)
end


nothing
