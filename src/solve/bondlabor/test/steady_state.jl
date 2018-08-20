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
    @test saved_ell ≈ m[:lstar].value
    @test saved_c   ≈ m[:cstar].value
    @test saved_η   ≈ m[:ηstar].value
    @test saved_μ   ≈ m[:μstar].value
    @test saved_β   ≈ m[:βstar].value
end
