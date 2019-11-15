using PDMats, Distributions
writing_output = false

####################################################################
# Testing Adaptive Φ Solution
####################################################################

file   = JLD2.jldopen("$path/../../reference/helpers_input.jld2", "r")
cloud  = old_to_new_cloud(DSGE.Cloud(read(file, "cloud")))
i_smc  = read(file, "i")
j_smc  = read(file, "j")
ϕ_prop = read(file, "phi_prop")
ϕ_n1   = read(file, "phi_n1")
proposed_fixed_schedule = read(file, "proposed_fixed_schedule")
tempering_target        = read(file, "tempering_target")
resampled_last_period   = read(file, "resampled_last_period")
close(file)

@everywhere Random.seed!(42)
test_ϕ_n, test_resampled_last_period, test_j, test_ϕ_prop = SMC.solve_adaptive_ϕ(cloud,
                                                                proposed_fixed_schedule,
                                                                i_smc, j_smc, ϕ_prop, ϕ_n1,
                                                                tempering_target,
                                                                resampled_last_period)
if writing_output
    jldopen("$path/../../reference/helpers_output.jld2", "w") do file
        write(file, "phi_n", test_ϕ_n)
        write(file, "resampled_last_period", test_resampled_last_period)
        write(file, "j", test_j)
        write(file, "phi_prop", test_ϕ_prop)
    end
    jldopen("$path/../../reference/helpers_output.jld2", "w") do file
        write(file, "phi_n", test_ϕ_n)
        write(file, "resampled_last_period", test_resampled_last_period)
        write(file, "j", test_j)
        write(file, "phi_prop", test_ϕ_prop)
    end
end

file = JLD2.jldopen("$path/../../reference/helpers_output.jld2", "r")
saved_ϕ_n = read(file, "phi_n")
saved_resampled_last_period = read(file, "resampled_last_period")
saved_j = read(file, "j")
saved_ϕ_prop = read(file, "phi_prop")
close(file)

####################################################################
@testset "Solve Adaptive Φ" begin
    @test test_ϕ_n ≈ saved_ϕ_n
    @test test_resampled_last_period == saved_resampled_last_period
    @test test_j == saved_j
    @test test_ϕ_prop == saved_ϕ_prop
end

####################################################################
# Testing MvNormal Mixture Draw Function
####################################################################
#=
file = JLD2.jldopen("reference/mutation_inputs.jld2")
saved_particles = read(file, "particles")
blocks_free     = read(file, "blocks_free")
blocks_all      = read(file, "blocks_all")
d               = read(file, "d")
α         = read(file, "α")
c         = read(file, "c")
close(file)
para = saved_particles[1].value
for (block_f, block_a) in zip(blocks_free, blocks_all)
    para_subset = para[block_a]
    d_subset    = MvNormal(d.μ[block_f], d.Σ.mat[block_f, block_f])
    test_θ_new, test_new_mix_density, test_old_mix_density = DSGE.mvnormal_mixture_draw(para_subset,
                                                                                   d_subset; c=c, α=α)
    JLD2.jldopen("reference/mvnormal_inputs.jld2", true, true, true, IOStream) do file
        write(file, "para_subset", para_subset)
        write(file, "d_subset", d_subset)
        write(file, "α", α)
        write(file, "c", c)
    end
    JLD2.jldopen("reference/mvnormal_output.jld2", true, true, true, IOStream) do file
        write(file, "θ_new", test_θ_new)
        write(file, "new_mix_density", test_new_mix_density)
        write(file, "old_mix_density", test_old_mix_density)
    end
end
=#

file = JLD2.jldopen("$path/../../reference/mvnormal_inputs.jld2")
para_subset = read(file, "para_subset")
d_subset    = read(file, "d_subset")
α           = read(file, "α")
c           = read(file, "c")
close(file)

test_θ_new, test_new_mix_density, test_old_mix_density = DSGE.mvnormal_mixture_draw(para_subset,
                                                                                   d_subset; c=c, α=α)

#=JLD2.jldopen("reference/mvnormal_output.jld2", true, true, true, IOStream) do file
    write(file, "θ_new", test_θ_new)
    write(file, "new_mix_density", test_new_mix_density)
    write(file, "old_mix_density", test_old_mix_density)
end=#

file = JLD2.jldopen("$path/../../reference/mvnormal_output.jld2")
    saved_θ_new = read(file, "θ_new")
    saved_new_mix_density = read(file, "new_mix_density")
    saved_old_mix_density = read(file, "old_mix_density")
close(file)

####################################################################
@testset "MvNormal Mixture Draw" begin
    @test test_θ_new == saved_θ_new
    @test test_new_mix_density == saved_new_mix_density
    @test test_old_mix_density == saved_old_mix_density
end

####################################################################
# Testing ESS Computation
####################################################################
#=
file = jldopen("reference/smc_sw_cloud_fix=true_blocks=3.jld2")
cloud = read(file, "cloud")
current_weights = read(file, "w")[:,3]
close(file)

n_part    = length(cloud.particles)
loglh     = [cloud.particles[i].loglh for i=1:n_part]
old_loglh = [cloud.particles[i].old_loglh for i=1:n_part]
ϕ_n       = 9.25022e-6
ϕ_n1      = 2.15769e-6

ess = DSGE.compute_ESS(loglh, current_weights, ϕ_n, ϕ_n1)

jldopen("reference/ess_inputs.jld2", true, true, true, IOStream) do file
    write(file, "loglh", loglh)
    write(file, "current_weights", current_weights)
    write(file, "ϕ_n", ϕ_n)
    write(file, "ϕ_n1", ϕ_n1)
    write(file, "old_loglh", old_loglh)
end

jldopen("reference/ess_output.jld2", true, true, true, IOStream) do file
    write(file, "ess", ess)
end
=#
file = JLD2.jldopen("$path/../../reference/ess_inputs.jld2")
loglh = read(file, "loglh")
current_weights = read(file, "current_weights")
ϕ_n = read(file, "ϕ_n")
ϕ_n1 = read(file, "ϕ_n1")
close(file)

file = JLD2.jldopen("$path/../../reference/ess_output.jld2")
saved_ESS = read(file, "ess")
close(file)

test_ESS = SMC.compute_ESS(loglh, current_weights, ϕ_n, ϕ_n1)
####################################################################
@testset "Compute ESS" begin
    @test test_ESS ≈ saved_ESS
end
