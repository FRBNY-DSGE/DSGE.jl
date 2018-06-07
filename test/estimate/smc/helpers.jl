# To be removed after running this test individually in the REPL successfully
using DSGE
using HDF5, JLD
import Base.Test: @test, @testset

file = jldopen("../../reference/helpers_input.jld", "r")
cloud = read(file, "cloud")
proposed_fixed_schedule = read(file, "proposed_fixed_schedule")
i = read(file, "i")
j = read(file, "j")
ϕ_prop = read(file, "phi_prop")
ϕ_n1 = read(file, "phi_n1")
tempering_target = read(file, "tempering_target")
resampled_last_period = read(file, "resampled_last_period")
close(file)

srand(42)
test_ϕ_n, test_resampled_last_period, test_j, test_ϕ_prop = DSGE.solve_adaptive_ϕ(cloud, proposed_fixed_schedule,
                                                                             i, j, ϕ_prop, ϕ_n1, tempering_target,
                                                                             resampled_last_period)

file = jldopen("../../reference/helpers_output.jld", "r")
saved_ϕ_n = read(file, "phi_n")
saved_resampled_last_period = read(file, "resampled_last_period")
saved_j = read(file, "j")
saved_ϕ_prop = read(file, "phi_prop")
close(file)

####################################################################
@testset "Solve Adaptive Φ" begin
    @test test_ϕ_n == saved_ϕ_n
    @test test_resampled_last_period == saved_resampled_last_period
    @test test_j == saved_j
    @test test_ϕ_prop == saved_ϕ_prop
end
