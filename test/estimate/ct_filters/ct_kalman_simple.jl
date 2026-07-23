using DSGE, Test, BenchmarkTools, LinearAlgebra, Random

# Test + benchmark for ct_kalman_simple (src/estimate/ct_filters/ct_kalman_simple.jl).


module CTSimpleScaffold
    using LinearAlgebra
    include(joinpath(@__DIR__, "..", "..", "..", "src", "estimate",
                     "ct_filters", "ct_kalman_simple.jl"))
end

Random.seed!(47)
T      = [-0.5 0.0; 0.0 -0.3]
Z      = [1.0 0.0]
Q      = 0.1 * Matrix{Float64}(I, 2, 2)
E      = reshape([0.05], 1, 1)
mean_0 = [0.0, 0.0]
var_0  = 0.1 * Matrix{Float64}(I, 2, 2)
n_data = 20
data_y = 0.1 * randn(n_data, 1)
dt     = fill(0.25, n_data)

out, run_err = nothing, nothing
try
    global out, run_err
    out = CTSimpleScaffold.ct_kalman_simple(T, Z, Q, E, mean_0, var_0, data_y, dt)
catch e
    global out, run_err
    run_err = e
end

@testset "ct_kalman_simple on a small stable system" begin
    @test run_err === nothing
    if run_err === nothing
        @test out isa Real
        @test isfinite(out)
    end
end

run_benchmarks = false
if run_benchmarks && run_err === nothing
    b = @benchmark CTSimpleScaffold.ct_kalman_simple($T, $Z, $Q, $E, $mean_0,
                                                     $var_0, $data_y, $dt)
    println("\nct_kalman_simple  time: ", BenchmarkTools.prettytime(median(b).time),
            "   memory: ", BenchmarkTools.prettymemory(median(b).memory))
end

nothing
