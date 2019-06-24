# To be removed after running this test individually in the REPL successfully
using DSGE
using HDF5, JLD2, Distributions, LinearAlgebra, PDMats
using Test

m = AnSchorfheide(testing = true)

@load "reference/smc.jld2" data

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

@load "reference/mutation_inputs.jld2" old_particles d blocks_free blocks_all ϕ_n ϕ_n1 c α old_data

function stack_values(particles::Vector{Particle}, field::Symbol)
    n_parts = length(particles)
    init_field_value = getfield(particles[1], field)
    stacked_values = Vector{typeof(init_field_value)}(undef, n_parts)

    for (i, particle) in enumerate(particles)
        stacked_values[i] = getfield(particle, field)
    end

    return stacked_values
end

Random.seed!(42)

new_particles = [mutation(m, data, old_particles[j], d, blocks_free, blocks_all, ϕ_n, ϕ_n1;
                 c = c, α = α, old_data = old_data) for j = 1:n_parts]

saved_particles = load("reference/mutation_outputs.jld2", "particles")

particle_fields = fieldnames(typeof(new_particles[1]))
@testset "Individual Particle Fields Post-Mutation" begin
    @test stack_values(new_particles, :weight) == stack_values(saved_particles, :weight)
    @test stack_values(new_particles, :keys) == stack_values(saved_particles, :keys)
    @test stack_values(new_particles, :value) == stack_values(saved_particles, :value)
    @test stack_values(new_particles, :loglh) ≈ stack_values(saved_particles, :loglh)
    @test stack_values(new_particles, :logpost) ≈ stack_values(saved_particles, :logpost)
    @test stack_values(new_particles, :old_loglh) == stack_values(saved_particles, :old_loglh)
    @test stack_values(new_particles, :accept) == stack_values(saved_particles, :accept)
end
