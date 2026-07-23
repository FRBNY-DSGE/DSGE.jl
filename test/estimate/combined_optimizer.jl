using DSGE, Test, BenchmarkTools, Random

# Tests + benchmark for combined_optimizer on a 2-D Rosenbrock (min 0 at [1,1]).
# (combined_optimizer.jl migrated to the Optim 2 API: LBFGS(), autodiff as an ADTypes optimize
#  kwarg, positional method + Optim.Options for SA, neighbor= keyword, round(...; digits=).)
rosenbrock(x::Vector) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
x0 = [-1.2, 1.0]

Random.seed!(47)
out = combined_optimizer(rosenbrock, copy(x0);
                         iterations = 100, max_cycles = 2, verbose = :none)

@testset "combined_optimizer on Rosenbrock" begin
    @test out.minimizer ≈ [1.0, 1.0] atol = 1e-2
    @test out.minimum   <  1e-4
end

# Benchmark (skipped while the function errors).
run_benchmarks = false
if run_benchmarks
    b = @benchmark combined_optimizer($rosenbrock, $(copy(x0));
                                      iterations = 100, max_cycles = 2,
                                      verbose = :none) evals = 1
    println("\ncombined_optimizer [Rosenbrock]  time: ",
            BenchmarkTools.prettytime(median(b).time),
            "   memory: ", BenchmarkTools.prettymemory(median(b).memory))
end

nothing
