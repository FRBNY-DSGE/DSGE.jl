using DSGE, Test, BenchmarkTools, LinearAlgebra, Random

# Test + benchmark for block_kalman_filter (loaded by DSGE).

Ttild  = [-0.5 0.0; 0.0 -0.3]
Rtild  = Matrix{Float64}(I, 2, 2)
Ctild  = [0.0, 0.0]
Qtild  = 0.1 * Matrix{Float64}(I, 2, 2)
Ztild  = [1.0 0.0]
Dtild  = [0.0]
Etild  = reshape([0.05], 1, 1)
M      = Matrix{Float64}(I, 2, 2)
Mtild  = Matrix{Float64}(I, 2, 2)
block_dims = [1, 0, 0, 1]
s_0tild = Vector{Float64}(undef, 0)
P_0tild = Matrix{Float64}(undef, 0, 0)
y = 0.1 * randn(1, 10)

Random.seed!(47)
out = block_kalman_filter(y, Ttild, Rtild, Ctild, Qtild, Ztild, Dtild, Etild,
                          M, Mtild, block_dims, s_0tild, P_0tild)

@testset "block_kalman_filter on a small 2-block system" begin
    loglh = out[1]
    @test length(loglh) == size(y, 2)
    @test all(isfinite, loglh)
end

run_benchmarks = false
if run_benchmarks
    b = @benchmark block_kalman_filter($y, $Ttild, $Rtild, $Ctild, $Qtild, $Ztild,
                                       $Dtild, $Etild, $M, $Mtild, $block_dims,
                                       $s_0tild, $P_0tild)
    println("\nblock_kalman_filter  time: ", BenchmarkTools.prettytime(median(b).time),
            "   memory: ", BenchmarkTools.prettymemory(median(b).memory))
end

nothing
