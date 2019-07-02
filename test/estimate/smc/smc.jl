using HDF5, DSGEModels, FileIO, JLD2, Test, Random, DSGE, DelimitedFiles
path = dirname(@__FILE__)

writing_output = false

m = AnSchorfheide()

save = normpath(joinpath(dirname(@__FILE__),"save"))
m <= Setting(:saveroot, save)

data = h5read("../../reference/smc.h5", "data")

m <= Setting(:n_particles, 400)
m <= Setting(:n_Φ, 100)
m <= Setting(:λ, 2.0)
m <= Setting(:n_smc_blocks, 1)
m <= Setting(:use_parallel_workers, true)
m <= Setting(:step_size_smc, 0.5)
m <= Setting(:n_mh_steps_smc, 1)
m <= Setting(:resampler_smc, :polyalgo)
m <= Setting(:target_accept, 0.25)

m <= Setting(:mixture_proportion, .9)
m <= Setting(:adaptive_tempering_target_smc, false)
m <= Setting(:resampling_threshold, .5)
m <= Setting(:smc_iteration, 0)
m <= Setting(:use_chand_recursion, true)

@everywhere Random.seed!(42)
smc(m, data, verbose = :none) # us.txt gives equiv to periods 95:174 in our current dataset

test_file = load(rawpath(m, "estimate", "smc_cloud.jld2"))
test_cloud  = test_file["cloud"]
test_w      = test_file["w"]
test_W      = test_file["W"]
test_z = test_file["z"]

if writing_output
    jldopen("../../reference/smc_cloud_fix=true.jld2", true, true, true, IOStream) do file
        write(file, "cloud", test_cloud)
        write(file, "w", test_w)
        write(file, "W", test_W)
        write(file, "z", test_z)
    end
end

saved_file = load("../../reference/smc_cloud_fix=true.jld2")
saved_cloud  = saved_file["cloud"]
saved_w      = saved_file["w"]
saved_W      = saved_file["W"]


####################################################################
cloud_fields = fieldnames(typeof(test_cloud))
@testset "ParticleCloud Fields: AnSchorf" begin
    @test @test_matrix_approx_eq DSGE.get_vals(test_cloud) DSGE.get_vals(saved_cloud)
    @test @test_matrix_approx_eq DSGE.get_loglh(test_cloud) DSGE.get_loglh(saved_cloud)
    @test length(test_cloud.particles) == length(saved_cloud.particles)
    @test test_cloud.tempering_schedule == saved_cloud.tempering_schedule
    @test test_cloud.ESS == saved_cloud.ESS
    @test test_cloud.stage_index == saved_cloud.stage_index
    @test test_cloud.n_Φ == saved_cloud.n_Φ
    @test test_cloud.resamples == saved_cloud.resamples
    @test test_cloud.c == saved_cloud.c
    @test test_cloud.accept == saved_cloud.accept
end

test_particle = test_cloud.particles[1]
saved_particle = saved_cloud.particles[1]
particle_fields = fieldnames(typeof(test_particle))
@testset "Individual Particle Fields Post-SMC: AnSchorf" begin
    @test test_particle.weight == saved_particle.weight
    @test test_particle.keys == saved_particle.keys
    @test test_particle.value == saved_particle.value
    @test test_particle.loglh == saved_particle.loglh
    @test test_particle.logpost == saved_particle.logpost
    @test test_particle.old_loglh == saved_particle.old_loglh
    @test test_particle.accept == saved_particle.accept
end

@testset "Weight Matrices: AnSchorf" begin
    @test @test_matrix_approx_eq test_w saved_w
    @test @test_matrix_approx_eq test_W saved_W
end
