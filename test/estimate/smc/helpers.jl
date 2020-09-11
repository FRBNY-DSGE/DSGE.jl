writing_output = false
if VERSION < v"1.5"
    ver = "111"
else 
    ver = "150"
end
@everywhere Random.seed!(42)

####################################################################
# Testing Adaptive Φ Solution
####################################################################
file   = JLD2.jldopen(joinpath(dirname(@__FILE__), "reference/solve_adaptive_phi.jld2"), "r")
cloud  = read(file, "cloud")
i_smc  = read(file, "i")
j_smc  = read(file, "j")
ϕ_prop = read(file, "phi_prop")
ϕ_n1   = read(file, "phi_n1")
proposed_fixed_schedule = read(file, "proposed_fixed_schedule")
tempering_target        = read(file, "tempering_target")
resampled_last_period   = read(file, "resampled_last_period")
close(file)

test_ϕ_n, test_resampled_last_period, test_j, test_ϕ_prop = SMC.solve_adaptive_ϕ(cloud,
                                                                proposed_fixed_schedule,
                                                                i_smc, j_smc, ϕ_prop, ϕ_n1,
                                                                tempering_target,
                                                                resampled_last_period)
if writing_output
    jldopen(joinpath(dirname(@__FILE__), "reference/helpers_output_version=" * ver * ".jld2"),
            true, true, true, IOStream) do file
        write(file, "phi_n", test_ϕ_n)
        write(file, "resampled_last_period", test_resampled_last_period)
        write(file, "j", test_j)
        write(file, "phi_prop", test_ϕ_prop)
    end
end

file = JLD2.jldopen(joinpath(dirname(@__FILE__), "reference/helpers_output_version=" 
                             * ver * ".jld2"), "r")
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
file = JLD2.jldopen(joinpath(dirname(@__FILE__),"reference/mvnormal_inputs.jld2"))
    para_subset = read(file, "para_subset")
    d_subset    = read(file, "d_subset")
    α           = read(file, "α")
    c           = read(file, "c")
close(file)

test_θ_new = SMC.mvnormal_mixture_draw(para_subset, d_subset; c=c, α=α)

if writing_output
    JLD2.jldopen(joinpath(dirname(@__FILE__),"reference/mvnormal_output_version=" * ver * ".jld2"), true, true, true, IOStream) do file
        write(file, "θ_new", test_θ_new)
    end
end

file = JLD2.jldopen(joinpath(dirname(@__FILE__),"reference/mvnormal_output_version=" * ver * ".jld2"))
    saved_θ_new = read(file, "θ_new")
close(file)

####################################################################
@testset "MvNormal Mixture Draw" begin
    @test test_θ_new == saved_θ_new
end


####################################################################
# Test: get_cov()
####################################################################
d = JLD2.jldopen(joinpath(dirname(@__FILE__),"reference/mutation_inputs.jld2"), "r") do file
    file["d"]
end
d_deg = DegenerateMvNormal(d.μ, d.Σ.mat)

####################################################################
@testset "Test: get_cov(d)" begin
    @test SMC.get_cov(d)     == d.Σ.mat
    @test SMC.get_cov(d_deg) == d_deg.σ
end


####################################################################
# Test: compute_proposal_densities()
####################################################################
file = JLD2.jldopen(joinpath(dirname(@__FILE__),"reference/proposal_densities_in.jld2"))
    para_draw   = read(file, "para_draw")
    para_subset = read(file, "para_subset")
    d_subset    = read(file, "d_subset")
    α           = read(file, "α")
    c           = read(file, "c")
close(file)

q0, q1 = SMC.compute_proposal_densities(para_draw, para_subset, d_subset; α = α,
                                        c = c)
if writing_output
    JLD2.jldopen(joinpath(dirname(@__FILE__),"reference/proposal_densities_output_version=" * ver * ".jld2"), true, true, true, IOStream) do file
        file["q0"] = q0
        file["q1"] = q1
    end
end

file = JLD2.jldopen(joinpath(dirname(@__FILE__),"reference/proposal_densities_output_version=" * ver * ".jld2"))
    saved_q0 = read(file, "q0")
    saved_q1 = read(file, "q1")
close(file)

####################################################################
@testset "Proposal densities" begin
    @test q0 ≈ saved_q0
    @test q1 ≈ saved_q1
end


####################################################################
# Testing ESS Computation
####################################################################
file = JLD2.jldopen(joinpath(dirname(@__FILE__),"reference/ess_inputs.jld2"))
    loglh           = read(file, "loglh")
    current_weights = read(file, "current_weights")
    ϕ_n             = read(file, "ϕ_n")
    ϕ_n1            = read(file, "ϕ_n1")
close(file)

test_ESS = SMC.compute_ESS(loglh, current_weights, ϕ_n, ϕ_n1)

if writing_output
    file = JLD2.jldopen(joinpath(dirname(@__FILE__),"reference/smc_sw_cloud_fix=true_blocks=3.jld2"), "r") 
    saved_cloud = read(file, "cloud")
    current_weights = read(file, "w")[:,3]
    close(file)

    n_part    = length(saved_cloud.particles)
    println(typeof(saved_cloud))
    loglh     = [saved_cloud.particles[i].loglh for i=1:n_part]
    old_loglh = [saved_cloud.particles[i].old_loglh for i=1:n_part]
    ϕ_n       = 9.25022e-6
    ϕ_n1      = 2.15769e-6

    JLD2.jldopen(joinpath(dirname(@__FILE__),"reference/ess_inputs.jld2"), true, true, true, IOStream) do file
        write(file, "loglh", loglh)
        write(file, "current_weights", current_weights)
        write(file, "ϕ_n", ϕ_n)
        write(file, "ϕ_n1", ϕ_n1)
        write(file, "old_loglh", old_loglh)
    end

    JLD2.jldopen(joinpath(dirname(@__FILE__),"reference/ess_output_version=" * ver * ".jld2"), true, true, true, IOStream) do file
        write(file, "ess", test_ESS)
    end
end

file = JLD2.jldopen(joinpath(dirname(@__FILE__),"reference/ess_output_version=" * ver * ".jld2"))
    saved_ESS = read(file, "ess")
close(file)

####################################################################
@testset "Compute ESS" begin
    @test test_ESS ≈ saved_ESS
end

####################################################################
# Testing Block Creation
####################################################################
m = AnSchorfheide()

free_para_inds = findall(x -> !x.fixed, m.parameters)
n_free_para    = length(free_para_inds)
n_blocks       = 3

test_blocks_free = SMC.generate_free_blocks(n_free_para, n_blocks)
test_blocks_all  = SMC.generate_all_blocks(test_blocks_free, free_para_inds)
test_blocks      = SMC.generate_param_blocks(length(m.parameters), n_blocks)

if writing_output
    JLD2.jldopen(joinpath(dirname(@__FILE__),"reference/helpers_blocking_version=" * ver * ".jld2"), true, true, true, IOStream) do file
        file["blocks_free"] = test_blocks_free
        file["blocks_all"]  = test_blocks_all
        file["blocks"]      = test_blocks
    end
end

saved_blocks_free = load(joinpath(dirname(@__FILE__),"reference/helpers_blocking_version=" * ver * ".jld2"), "blocks_free")
saved_blocks_all  = load(joinpath(dirname(@__FILE__),"reference/helpers_blocking_version=" * ver * ".jld2"), "blocks_all")
saved_blocks      = load(joinpath(dirname(@__FILE__),"reference/helpers_blocking_version=" * ver * ".jld2"), "blocks")

####################################################################
@testset "Mutation block generation" begin
    @test test_blocks_free == saved_blocks_free
    @test test_blocks_all  == saved_blocks_all
    @test test_blocks      == saved_blocks
end
