using Base.Test
using JLD

path = dirname(@__FILE__)

file = jldopen("$path/reference/solve.jld", "r")
saved_gx2  = read(file, "gx2")
saved_hx2  = read(file, "hx2")
saved_gx   = read(file, "gx")
saved_hx   = read(file, "hx")
close(file)

@testset "Check solve outputs" begin
    @test saved_gx2 ≈ gx2
    @test saved_hx2 ≈ hx2
    @test saved_gx  ≈ gx
    @test saved_hx  ≈ hx
end
