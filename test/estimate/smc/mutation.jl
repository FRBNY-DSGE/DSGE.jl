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
m <= Setting(:resampler_smc, :multinomial)
m <= Setting(:target_accept, 0.25)

m <= Setting(:mixture_proportion, .9)
m <= Setting(:tempering_target, 0.95)
m <= Setting(:resampling_threshold, .5)
m <= Setting(:use_fixed_schedule, true)

n_parts = get_setting(m, :n_particles)

file = jldopen("../../reference/mutation_inputs.jld", "r")
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

function stack_values(particles::Vector{Particle}, field::Symbol)
    n_parts = length(particles)
    init_field_value = getfield(particles[1], field)
    stacked_values = Vector{typeof(init_field_value)}(n_parts)

    for (i, particle) in enumerate(particles)
        stacked_values[i] = getfield(particle, field)
    end

    return stacked_values
end

srand(42)

new_particles = [mutation(m, data, old_particles[j], d, blocks_free, blocks_all, ϕ_n, ϕ_n1;
                 c = c, α = α, old_data = old_data) for j = 1:n_parts]

saved_particles = load("../../reference/mutation_outputs.jld", "particles")

particle_fields = fieldnames(new_particles[1])
@testset "Individual Particle Fields Post-Mutation" begin
    @test stack_values(new_particles, :weight) == stack_values(saved_particles, :weight)
    for field in particle_fields
        new_particles_field = stack_values(new_particles, field)
        saved_particles_field = stack_values(saved_particles, field)

        @test new_particles_field == saved_particles_field
    end
end
