using DSGE, Test, BenchmarkTools, LinearAlgebra, Random

# Test + benchmark for nearest_spd(A) in src/estimate/nearest_spd.jl
# (Higham's nearest symmetric positive semidefinite matrix).
#
# The source calls chol(Ahat), but `chol` was removed in Julia 1.0 (-> cholesky).
# As written, chol throws UndefVarError on every iteration of the PD-check loop;
# the bare catch swallows it, so chol_success never flips -> the loop never exits
# (infinite hang, not a clean error). Shim `chol` into DSGE so it resolves — the
# real fix is to replace chol(Ahat) with cholesky(Ahat) in the source.
Core.eval(DSGE, :(chol(A) = LinearAlgebra.cholesky(A)))

Random.seed!(47)

@testset "nearest_spd: already-SPD matrix returned ~unchanged" begin
    M = randn(5, 5)
    A = M' * M + 5I                 # symmetric positive definite
    Ahat = nearest_spd(Matrix(A))
    @test Ahat ≈ A
end

@testset "nearest_spd: arbitrary matrix -> symmetric PD" begin
    A = randn(6, 6)                 # non-symmetric, generally indefinite
    Ahat = nearest_spd(A)
    @test Ahat ≈ Ahat'             # symmetric
    @test isposdef(Ahat)           # loop perturbs until cholesky succeeds
end

@testset "nearest_spd: non-square throws" begin
    @test_throws ErrorException nearest_spd(randn(3, 4))
end

run_benchmarks = false
if run_benchmarks
    A = randn(20, 20)
    b = @benchmark nearest_spd($A)
    println("\nnearest_spd (20x20)  time: ", BenchmarkTools.prettytime(median(b).time),
            "   memory: ", BenchmarkTools.prettymemory(median(b).memory))
end

nothing
