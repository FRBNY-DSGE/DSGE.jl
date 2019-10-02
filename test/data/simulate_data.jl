m = AnSchorfheide()
system = compute_system(m)

Random.seed!(1793)
s₀, P₀ = DSGE.init_stationary_states(system[:TTT], system[:RRR], system[:CCC], system[:QQ])
s = Matrix{Float64}(undef, size(system[:TTT],1), 2)
ϵ = rand(DegenerateMvNormal(zeros(size(system[:RRR],2)), sqrt.(system[:QQ])), 2)
s[:,1] = s₀
s[:,2] = system[:CCC] + system[:TTT] * s₀ + system[:RRR] * ϵ[:,2]
Random.seed!(1793)
sout, ϵout = simulate_states(system[:TTT], system[:RRR], system[:CCC],
                system[:QQ], burnin = 1,
                n_periods = 1)
sout1, ϵout1 = simulate_states(system[:TTT], system[:RRR], system[:CCC],
                system[:QQ], burnin = 1,
                n_periods = 1, ϵ = ϵ)
@testset "Check data can be properly simulated" begin
    @test @test_matrix_approx_eq sout s
    @test @test_matrix_approx_eq ϵout ϵ
    @test @test_matrix_approx_eq sout1 s
    @test @test_matrix_approx_eq ϵout1 ϵ
end
