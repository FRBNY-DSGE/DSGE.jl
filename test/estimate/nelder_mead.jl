using DSGE, Test, BenchmarkTools, Random
import Optim

# Tests + benchmark for src/estimate/nelder_mead.jl.
#   - nelder_mead(fcn, x0): uses Optim's positional method + Options API; we check the
#     objective improves substantially on Rosenbrock.
#   - MatlabSimplexer constructors.
#   - Optim.simplexer(MatlabSimplexer, x): n+1 distinct vertices (copies the initial point;
#     perturbs each coordinate, matching MATLAB fminsearch).

#-----------------------------------------------------------------
# nelder_mead
#-----------------------------------------------------------------
rosenbrock(x::Vector) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
x0 = [-1.2, 1.0]

Random.seed!(47)
out = DSGE.nelder_mead(rosenbrock, copy(x0); iterations = 1000)

@testset "nelder_mead on Rosenbrock" begin
    @test isfinite(out.minimum)
    @test out.minimum < 1.0           # big improvement from f(x0) ≈ 24.2
end

#-----------------------------------------------------------------
# MatlabSimplexer constructors
#-----------------------------------------------------------------
@testset "MatlabSimplexer constructors" begin
    s = DSGE.MatlabSimplexer()
    @test s.a == 0.00025
    @test s.b == 0.05

    s2 = DSGE.MatlabSimplexer(1.0, 2.0)
    @test s2.a == 1.0
    @test s2.b == 2.0

    s3 = DSGE.MatlabSimplexer(a = 0.1, b = 0.2)
    @test s3.a == 0.1
    @test s3.b == 0.2
end

#-----------------------------------------------------------------
# Optim.simplexer(MatlabSimplexer, x)
#-----------------------------------------------------------------
@testset "Optim.simplexer(MatlabSimplexer, x)" begin
    x = [1.0, 2.0, 3.0]
    simplex = Optim.simplexer(DSGE.MatlabSimplexer(), x)
    @test length(simplex) == length(x) + 1          # n+1 vertices
    @test allunique(simplex)                         # each vertex is a distinct, perturbed copy
end

################
# Benchmarking #
################
run_benchmarks = false
if run_benchmarks
    results = Tuple{String, Any}[]

    b_nm = @benchmark DSGE.nelder_mead($rosenbrock, $(copy(x0)); iterations = 1000)
    push!(results, ("nelder_mead", b_nm))

    b_simplex = @benchmark Optim.simplexer(DSGE.MatlabSimplexer(), x) setup = (x = [1.0, 2.0, 3.0, 4.0])
    push!(results, ("Optim.simplexer (4-vec)", b_simplex))

    println("\n===== estimate/nelder_mead benchmark results =====")
    for (name, b) in results
        println(rpad(name, 26), " time: ", BenchmarkTools.prettytime(median(b).time),
                "   memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end

nothing
