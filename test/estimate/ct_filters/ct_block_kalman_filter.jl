using DSGE, Test, BenchmarkTools, LinearAlgebra, Random

# Test + benchmark for ct_block_kalman_filter (NOT loaded by DSGE; commented out in
# DSGE.jl). Included directly here in an isolated module.
# Smoke/regression test: finite per-period log-likelihoods on a small stable system.

module CTBlockScaffold
    import DSGE: init_stationary_states
    using LinearAlgebra
    include(joinpath(@__DIR__, "..", "..", "..", "src", "estimate",
                     "ct_filters", "ct_block_kalman_filter.jl"))
end

T = [-0.5 0.0; 0.0 -0.3]
R = Matrix{Float64}(I, 2, 2)
C = [0.0, 0.0]
Q = 0.1 * Matrix{Float64}(I, 2, 2)
Z = [1.0 0.0]
D = [0.0]
E = reshape([0.05], 1, 1)
s_0 = zeros(2)
P_0 = Matrix{Float64}(I, 2, 2)
y = 0.1 * randn(1, 10)

Random.seed!(47)
out = CTBlockScaffold.ct_block_kalman_filter(y, T, R, C, Q, Z, D, E;
                                             n_simulate_states = 1, s_0 = s_0, P_0 = P_0)

@testset "ct_block_kalman_filter on a small system" begin
    loglh = out[1]
    @test length(loglh) == size(y, 2)
    @test all(isfinite, loglh)
end

run_benchmarks = false
if run_benchmarks
    b = @benchmark CTBlockScaffold.ct_block_kalman_filter($y, $T, $R, $C, $Q, $Z, $D, $E;
                                                          n_simulate_states = 1,
                                                          s_0 = $s_0, P_0 = $P_0)
    println("\nct_block_kalman_filter  time: ", BenchmarkTools.prettytime(median(b).time),
            "   memory: ", BenchmarkTools.prettymemory(median(b).memory))
end

nothing
