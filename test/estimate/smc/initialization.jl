# To be removed after running this test individually in the REPL successfully
using DSGE
using HDF5, JLD
import Base.Test: @test, @testset

path = dirname(@__FILE__)

m = AnSchorfheide()

saveroot = normpath(joinpath(dirname(@__FILE__),"..","..","save"))
m <= Setting(:saveroot, saveroot)

data = h5read("../../reference/smc.h5", "data")

m <= Setting(:n_particles, 400)
m <= Setting(:n_Φ, 100)
m <= Setting(:λ, 2.0)
m <= Setting(:n_smc_blocks, 1)
m <= Setting(:use_parallel_workers, false)
m <= Setting(:step_size_smc, 0.5)
m <= Setting(:n_mh_steps_smc, 1)
m <= Setting(:resampler_smc, :polyalgo)
m <= Setting(:target_accept, 0.25)

m <= Setting(:mixture_proportion, .9)
m <= Setting(:tempering_target, 0.95)
m <= Setting(:resampling_threshold, .5)
m <= Setting(:use_fixed_schedule, true)

####################################################################
test_init_cloud = ParticleCloud(m, get_setting(m, :n_particles))
srand(42)
DSGE.initial_draw!(m, data, test_init_cloud)
saved_init_cloud = load("../../reference/initial_draw.jld", "cloud")
@testset "Initial draw" begin
    @test @test_matrix_approx_eq DSGE.get_vals(test_init_cloud) DSGE.get_vals(saved_init_cloud)
    @test @test_matrix_approx_eq DSGE.get_loglh(test_init_cloud) DSGE.get_loglh(saved_init_cloud)
end

