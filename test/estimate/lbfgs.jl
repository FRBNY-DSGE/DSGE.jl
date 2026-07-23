using DSGE, Test, BenchmarkTools, Random

# lbfgs uses a bare LBFGS() that isn't in DSGE's namespace (DSGE only imports
# `optimize` from Optim). Bring LBFGS into DSGE so the function resolves — this is
# the same one-line import the fix would add to DSGE; remove once that's in.
Core.eval(DSGE, :(import Optim: LBFGS))

rosenbrock(x::Vector) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
x0 = [-1.2, 1.0]

Random.seed!(47)
out = lbfgs(rosenbrock, copy(x0); iterations = 1000)

@testset "lbfgs on Rosenbrock" begin
    @test out.minimizer ≈ [1.0, 1.0] atol = 1e-3
    @test out.minimum   <  1e-6
end

run_benchmarks = false
if run_benchmarks
    b = @benchmark lbfgs($rosenbrock, $(copy(x0)); iterations = 1000)
    println("\nlbfgs [Rosenbrock]  time: ", BenchmarkTools.prettytime(median(b).time),
            "   memory: ", BenchmarkTools.prettymemory(median(b).memory))
end

nothing
