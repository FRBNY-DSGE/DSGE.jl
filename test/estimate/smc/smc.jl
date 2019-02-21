# To be removed after running this test individually in the REPL successfully
using DSGE
using JLD2
using Test

m = AnSchorfheide()

m <= Setting(:saveroot, tempdir())

@load "reference/smc.jld2" data

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

Random.seed!(42)

smc(m, data, verbose = :none) # us.txt gives equiv to periods 95:174 in our current dataset

file = JLD2.jldopen(rawpath(m, "estimate", "smc_cloud.jld2"), "r")
test_cloud = read(file, "cloud")
test_w     = read(file, "w")
test_W     = read(file, "W")
close(file)

file = JLD2.jldopen("reference/smc_test_cloud.jld2", "r")
saved_cloud = read(file, "cloud")
saved_w     = read(file, "w")
saved_W     = read(file, "W")
close(file)

####################################################################
cloud_fields = fieldnames(typeof(test_cloud))
@testset "ParticleCloud Fields" begin
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
@testset "Individual Particle Fields Post-SMC" begin
    @test test_particle.weight == saved_particle.weight
    @test test_particle.keys == saved_particle.keys
    @test test_particle.value == saved_particle.value
    @test test_particle.loglh == saved_particle.loglh
    @test test_particle.logpost == saved_particle.logpost
    @test test_particle.old_loglh == saved_particle.old_loglh
    @test test_particle.accept == saved_particle.accept
end

@testset "Weight Matrices" begin
    @test @test_matrix_approx_eq test_w saved_w
    @test @test_matrix_approx_eq test_W saved_W
end
