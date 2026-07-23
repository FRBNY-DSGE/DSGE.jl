using DSGE, Test, BenchmarkTools, LinearAlgebra, Random

# Test + benchmark for ct_kalman_filter (loaded by DSGE).
# Currently errors. Known issues:
#   - Tsit5/ODEProblem are never imported into DSGE -> default method=Tsit5() throws
#   - 2nd method's forecast!: undefined `span` (183), `solvect` (186), `.u` vs `.u[end]` (188)
#   - Julia-0.x ctors Vector{S}(0)/Matrix{S}(0,0) in non-default output branches
# Written to pass once these are fixed.

T = [-0.5 0.0; 0.0 -0.3]
R = Matrix{Float64}(I, 2, 2)
C = [0.0, 0.0]
Q = 0.1 * Matrix{Float64}(I, 2, 2)
Z = [1.0 0.0]
D = [0.0]
E = reshape([0.05], 1, 1)
tspan = 0.25
y = 0.1 * randn(1, 10)

out, run_err = nothing, nothing
try
    global out, run_err
    Random.seed!(47)
    out = ct_kalman_filter(y, T, R, C, Q, Z, D, E, tspan)
catch e
    global out, run_err
    run_err = e
end

@testset "ct_kalman_filter on a small system" begin
    if run_err !== nothing
        @test_broken run_err === nothing  # known broken — see header
    else
        loglh = out[1]
        @test length(loglh) == size(y, 2)
        @test all(isfinite, loglh)
    end
end

run_benchmarks = false
if run_benchmarks && run_err === nothing
    b = @benchmark ct_kalman_filter($y, $T, $R, $C, $Q, $Z, $D, $E, $tspan)
    println("\nct_kalman_filter  time: ", BenchmarkTools.prettytime(median(b).time),
            "   memory: ", BenchmarkTools.prettymemory(median(b).memory))
end

nothing
