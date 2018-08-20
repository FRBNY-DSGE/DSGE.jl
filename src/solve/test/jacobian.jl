using Base.Test
using JLD

path = dirname(@__FILE__)

file = jldopen("$path/reference/jacobian.jld", "r")
saved_JJ  = read(file, "JJ")
close(file)

@testset "Check jacobian outputs" begin
    @test saved_JJ ≈ JJ
end
