using HDF5
using Test
using DSGE

m = PoolModel()
Φ1, F_ϵ1, F_λ1 = transition(m)
Φ2, F_ϵ2, F_λ2 = solve(m)

@testset "Check state-space system matches reference" begin
    @test Φ1 == Φ2
    @test F_ϵ1 == F_ϵ2
    @test F_λ1 == F_λ2
end
