m = AnSchorfheide()
system = compute_system(m)
Φ, Ψ, F_ϵ, F_u = DSGE.compute_system_function(system)
zero_sys = DSGE.zero_system_constants(system)

@testset "Check System, Transition, Measurement, and PseduoMeasurement access functions" begin
    @test typeof(system) == System{Float64}
    @test typeof(system[:transition]) == Transition{Float64}
    @test typeof(system[:measurement]) == Measurement{Float64}
    @test typeof(system[:pseudo_measurement]) == PseudoMeasurement{Float64}
    @test typeof(system[:transition][:TTT]) == Matrix{Float64}
    @test typeof(system[:transition][:RRR]) == Matrix{Float64}
    @test typeof(system[:transition][:CCC]) == Vector{Float64}
    @test typeof(system[:measurement][:ZZ]) == Matrix{Float64}
    @test typeof(system[:measurement][:DD]) == Vector{Float64}
    @test typeof(system[:measurement][:QQ]) == Matrix{Float64}
    @test typeof(system[:measurement][:EE]) == Matrix{Float64}
    @test typeof(system[:pseudo_measurement][:ZZ_pseudo]) == Matrix{Float64}
    @test typeof(system[:pseudo_measurement][:DD_pseudo]) == Vector{Float64}
end

@testset "Check miscellaneous functions acting on System types" begin
    @test sum(zero_sys[:CCC]) == 0.
    @test sum(zero_sys[:DD]) == 0.
    @test sum(zero_sys[:DD_pseudo]) == 0.
    @test Φ(ones(size(system[:transition][:TTT],2)), zeros(size(system[:transition][:RRR],2))) ==
        system[:transition][:TTT] * ones(size(system[:transition][:TTT],2))
    @test Ψ(ones(size(system[:transition][:TTT],2))) ==
        system[:measurement][:ZZ] * ones(size(system[:transition][:TTT],2)) + system[:measurement][:DD]
    @test F_ϵ.μ == zeros(3)
    @test F_ϵ.Σ.mat == system[:measurement][:QQ]
    @test F_u.μ == zeros(3)
    @test F_u.Σ.mat == system[:measurement][:EE]
end

nothing
