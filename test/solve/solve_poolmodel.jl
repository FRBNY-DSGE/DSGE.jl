using HDF5
using Test
using DSGE
using BenchmarkTools

m = PoolModel()
Φ1, F_ϵ1, F_λ1 = transition(m)
Φ2, F_ϵ2, F_λ2 = solve(m)

run_benchmarks = false

if run_benchmarks
    b_solve = @benchmark solve($m)

    println("\n===== solve (PoolModel) benchmark results =====")
    println(rpad("solve", 18), " time: ", rpad(BenchmarkTools.prettytime(median(b_solve).time), 12),
            "memory: ", BenchmarkTools.prettymemory(median(b_solve).memory))
end

@testset "Check state-space system matches reference" begin
    @test Φ1 == Φ2
    @test F_ϵ1 == F_ϵ2
    @test F_λ1 == F_λ2
end
