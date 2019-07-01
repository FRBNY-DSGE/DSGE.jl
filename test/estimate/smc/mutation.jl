# To be removed after running this test individually in the REPL successfully
using DSGE
using HDF5, JLD, JLD2, Random, Distributions, PDMats
import Test: @test, @testset

write_test_output = false

path = dirname(@__FILE__)

m = AnSchorfheide()

save = normpath(joinpath(dirname(@__FILE__),"save"))
m <= Setting(:saveroot, saveroot)

data = h5read("../../reference/smc.h5", "data")

m <= Setting(:n_particles, 400)
m <= Setting(:n_Φ, 100)
m <= Setting(:λ, 2.0)
m <= Setting(:n_smc_blocks, 1)
m <= Setting(:use_parallel_workers, false)
m <= Setting(:step_size_smc, 0.5)
m <= Setting(:n_mh_steps_smc, 1)
m <= Setting(:resampler_smc, :multinomial)
m <= Setting(:target_accept, 0.25)

m <= Setting(:mixture_proportion, .9)
m <= Setting(:tempering_target, 0.95)
m <= Setting(:resampling_threshold, .5)
m <= Setting(:use_fixed_schedule, true)

n_parts = get_setting(m, :n_particles)

file = JLD2.jldopen("../../reference/mutation_inputs.jld2", "r")
old_particles = read(file, "particles")
d = read(file, "d")
blocks_free = read(file, "blocks_free")
blocks_all = read(file, "blocks_all")
ϕ_n = read(file, "ϕ_n")
ϕ_n1 = read(file, "ϕ_n1")
c = read(file, "c")
α = read(file, "α")
old_data = read(file, "old_data")
close(file)

old_part_cloud = DSGE.vector_particles_to_cloud(m, old_particles)
Random.seed!(42)
new_particles = [mutation(m, data, old_part_cloud.particles[j, :], d.μ, Matrix(d.Σ),
                          blocks_free, blocks_all, ϕ_n, ϕ_n1;
                          c = c, α = α, old_data = old_data) for j = 1:n_parts]

if write_test_output
    #=JLD.jldopen("reference/mutation_outputs.jld", "w") do file
    write(file, "particles", new_particles)
    end =#
    JLD2.jldopen("../../reference/mutation_outputs.jld2", "w") do file
        write(file, "particles", new_particles)
    end
end

saved_particles = load("../../reference/mutation_outputs.jld2", "particles")

@testset "Test mutation outputs, particle by particle" begin
    for i = 1:length(saved_particles)
        @test isapprox(saved_particles[i], new_particles[i], nans = true)
    end
end
