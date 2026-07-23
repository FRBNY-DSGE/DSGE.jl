using DSGE, Test, BenchmarkTools, LinearAlgebra, Random

# Scaffold for filter(m::AbstractCTModel, data, system, ...) in src/estimate/filter_hank.jl.
# Broken + orphaned (never called anywhere; the CT/HANK model tests are disabled).
# The chain breaks before filter even runs; known issues, in order encountered:
#   - KrusellSmithCT() ctor: default custom_settings is Vector{Setting{Bool}}, but the
#     inner ctor expects Vector{Setting} (type invariance) -> MethodError. CT models
#     can't be constructed, which is why their tests are disabled.
#   - filter_hank.jl:47 eye(TTTs): `eye` is eye(::Integer) only -> MethodError on a Matrix
#   - filter_hank.jl:39 defaults Vector{S}(0)/Matrix{S}(0,0): Julia-0.x ctors
# Everything is wrapped so the first failure is captured. @test_broken until fixed.

# Cheap dummy System (does not depend on the model).
Ns, Ne, Ny = 3, 1, 1
trans  = DSGE.Transition(randn(Ns, Ns), randn(Ns, Ne))
meas   = DSGE.Measurement(randn(Ny, Ns), randn(Ny),
                          0.1 * Matrix{Float64}(I, Ne, Ne),
                          0.05 * Matrix{Float64}(I, Ny, Ny))
system = DSGE.System(trans, meas)
data   = randn(Ny, 10)
s_0    = Vector{Float64}(undef, 0)
P_0    = Matrix{Float64}(undef, 0, 0)

m, out, run_err = nothing, nothing, nothing
try
    global m, out, run_err
    Random.seed!(47)
    m   = KrusellSmithCT()
    out = DSGE.filter(m, data, system, s_0, P_0)
catch e
    global m, out, run_err
    run_err = e
end

@testset "filter(::AbstractCTModel) on a small dummy system" begin
    if run_err !== nothing
        @test_broken run_err === nothing  # known broken — see header
    else
        @test out isa Kalman
    end
end

run_benchmarks = false
if run_benchmarks && run_err === nothing
    b = @benchmark DSGE.filter($m, $data, $system, $s_0, $P_0)
    println("\nfilter(::AbstractCTModel)  time: ", BenchmarkTools.prettytime(median(b).time),
            "   memory: ", BenchmarkTools.prettymemory(median(b).memory))
end

nothing
