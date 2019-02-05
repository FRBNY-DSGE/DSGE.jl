# To be removed after running this test individually in the REPL successfully
using DSGE
using HDF5, JLD, JLD2
import Base.Test: @test, @testset

writing_output = false

####################################################################
# Testing Adaptive Φ Solution
####################################################################

file   = JLD2.jldopen("reference/helpers_input.jld2", "r")
cloud  = read(file, "cloud")
i      = read(file, "i")
j      = read(file, "j")
ϕ_prop = read(file, "phi_prop")
ϕ_n1   = read(file, "phi_n1")
proposed_fixed_schedule = read(file, "proposed_fixed_schedule")
tempering_target        = read(file, "tempering_target")
resampled_last_period   = read(file, "resampled_last_period")
close(file)

@everywhere srand(42)
test_ϕ_n, test_resampled_last_period, test_j, test_ϕ_prop = DSGE.solve_adaptive_ϕ(cloud, proposed_fixed_schedule,
                                                                             i, j, ϕ_prop, ϕ_n1, tempering_target,
                                                                             resampled_last_period)
if writing_output
    JLD2.jldopen("reference/helpers_output.jld2", "w") do file
        write(file, "phi_n", test_ϕ_n)
        write(file, "resampled_last_period", test_resampled_last_period)
        write(file, "j", test_j)
        write(file, "phi_prop", test_ϕ_prop)
    end
    JLD2.jldopen("reference/helpers_output.jld2", "w") do file
        write(file, "phi_n", test_ϕ_n)
        write(file, "resampled_last_period", test_resampled_last_period)
        write(file, "j", test_j)
        write(file, "phi_prop", test_ϕ_prop)
    end
end

file = JLD2.jldopen("reference/helpers_output.jld2", "r")
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
JLD2.jldopen("reference/mutation_inputs.jld2") do file
    saved_particles = read(file, "particles")
    saved_α         = read(file, "α")
    saved_c         = read(file, "c")
end
DSGE.solve_adaptive_ϕ(cloud, proposed_fixed_schedule,
                      i, j, ϕ_prop, ϕ_n1, tempering_target,
                      resampled_last_period)

####################################################################
@testset "MvNormal Mixture Draw" begin
    @test test_θ_new == saved_θ_new
    @test test_new_mix_density == saved_new_mix_density
    @test test_old_mix_density == saved_old_mix_density
end


####################################################################
# Testing ESS Computation
####################################################################
JLD2.jldopen("reference/smc_sw_cloud_fix=true_blocks=3.jld2") do file
end

compute_ESS(loglh, current_weights, ϕ_n, ϕ_n1, old_loglh)


####################################################################
@testset "Compute ESS" begin
    @test test_ESS == saved_ESS
end
