using DSGE
using HDF5
using BenchmarkTools

path = dirname(@__FILE__)

# Test in generic case
# Global minimum at (a, a^2)
a = 1
b = 100
function rosenbrock_csminwel(x::Vector)
    return (a-x[1])^2.0 + b*(x[2]-x[1]^2.0)^2.0
end

# this is the actual gradient of the Rosenbrock function
function rosenbrock_grad(x::Array)
    dr = similar(x)
    dr[1] = -2*(a-x[1]) - 2b*(2x[1])*(x[2]-x[1]^2.0)
    dr[2] = 2b*(x[2]-x[1]^2.0)
    badg = false
    return dr, badg
end

# A really bad guess
x_init = [10.0, -9.0]
x_expected = [a, a^2.0]

@testset "Check csminwel gradient calculation for correctness" begin
    res_real_grad,  = csminwel(rosenbrock_csminwel, rosenbrock_grad, x_init)
    @test x_expected ≈ res_real_grad.minimizer
    res_numeric_grad,  = csminwel(rosenbrock_csminwel, x_init)
    @test x_expected ≈ res_numeric_grad.minimizer atol=1.0e-8
end

################
# Benchmarking #
################
# Flip to true to run; off by default. Pure optimization, no FRED API.
run_benchmarks = false

if run_benchmarks
    b_analytic = @benchmark csminwel(rosenbrock_csminwel, rosenbrock_grad, $x_init)
    b_numeric  = @benchmark csminwel(rosenbrock_csminwel, $x_init)

    println("\n===== estimate/csminwel benchmark results =====")
    for (name, b) in [("analytic gradient", b_analytic),
                      ("numeric gradient ", b_numeric)]
        println(name, "  time:   ", BenchmarkTools.prettytime(median(b).time),
                "   memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end

nothing
