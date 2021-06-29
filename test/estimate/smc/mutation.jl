write_test_output = false
if VERSION < v"1.5"
    ver = "111"
elseif VERSION < v"1.6"
    ver = "150"
else
    ver = "160"
end

path = dirname(@__FILE__)

m = AnSchorfheide()

save = normpath(joinpath(dirname(@__FILE__),"save"))
m <= Setting(:saveroot, saveroot)

data = h5read(joinpath(path, "reference/smc.h5"), "data")

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

file = JLD2.jldopen(joinpath(path, "reference/mutation_inputs.jld2"), "r")
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
function my_likelihood(parameters::ParameterVector, data::Matrix{Float64})::Float64
    DSGE.update!(m, parameters)
    DSGE.likelihood(m, data; sampler = false, catch_errors = true,
                    use_chand_recursion = true, verbose = :low)
end

Random.seed!(42)
new_particles = [SMC.mutation(my_likelihood, m.parameters, data,
                              old_part_cloud.particles[j, vcat(1:16, 18:22)], # somehow the parameter vector got extended to 22, with index 17 being an undefined Float64
                              d.μ, Matrix(d.Σ),
                              16, # 16 because this is what saved test output had (but really AS has only 13 free params)
                              blocks_free, blocks_all, ϕ_n, ϕ_n1;
                              c = c, α = α, old_data = old_data) for j = 1:n_parts]

if write_test_output
    JLD2.jldopen(joinpath(path, "reference/mutation_outputs_version=" * ver * ".jld2"), "w") do file
        write(file, "particles", new_particles)
    end
end

saved_particles = load(joinpath(path, "reference/mutation_outputs_version=" * ver * ".jld2"), "particles")

@testset "Test mutation outputs, particle by particle" begin
    for i = 1:length(saved_particles)
        @test isapprox(saved_particles[i], new_particles[i], nans = true)
    end
end
