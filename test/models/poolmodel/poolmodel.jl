###########################################################################
# Set up for testing PoolModel instantiation
###########################################################################
# filepath = pwd()
filepath = dirname(@__FILE__)

@testset "Check constructors with dynamic, equal, and static weight options" begin
    pm = PoolModel()
    @test typeof(pm) == PoolModel{Float64}
    pm = PoolModel(weight_type = :dynamic)
    @test typeof(pm) == PoolModel{Float64}
    @test get_setting(pm, :weight_type) == :dynamic
    pm = PoolModel(weight_type = :equal)
    @test typeof(pm) == PoolModel{Float64}
    @test get_setting(pm, :weight_type) == :equal
    pm = PoolModel(weight_type = :static)
    @test typeof(pm) == PoolModel{Float64}
    @test get_setting(pm, :weight_type) == :static
end

# Check compute_system, transition and measurement functions apply to PoolModel
@testset "Check solve and statespace functions apply to PoolModel" begin
    pm = PoolModel()
    Φ1, Ψ1, F_ϵ1, F_u1, F_λ1 = compute_system(pm)
    Φ2, F_ϵ2, F_λ2 = transition(pm)
    Ψ2, F_u2 = measurement(pm)

    @test Φ1 == Φ2
    @test Ψ1 == Ψ2
    @test F_ϵ1 == F_ϵ2
    @test F_λ1 == F_λ2
    @test F_u1 == F_u2
end

nothing
