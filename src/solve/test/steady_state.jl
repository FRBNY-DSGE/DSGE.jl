using Base.Test
using JLD

path = dirname(@__FILE__)

file = jldopen("$path/reference/steady_state.jld", "r")
saved_ell  = read(file, "ell")
saved_c    = read(file, "c")
saved_η    = read(file, "eta")
saved_μ    = read(file, "mu")
saved_β    = read(file, "beta")
close(file)

@testset "Check steady state outputs" begin
    @test saved_ell ≈ ell
    @test saved_c   ≈ c
    @test saved_η   ≈ η
    @test saved_μ   ≈ μ
    @test saved_β   ≈ β
end
