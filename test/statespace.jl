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


@testset "VAR approximation of state space" begin
    m = Model1002("ss10"; custom_settings =
                  Dict{Symbol,Setting}(:add_laborshare_measurement =>
                                       Setting(:add_laborshare_measurement, true)))
    system = compute_system(m)
    system = compute_system(m, system; observables = [:obs_hours, :obs_gdpdeflator,
                                                      :laborshare_t, :NominalWageGrowth],
                            shocks = collect(keys(m.exogenous_shocks)))
    yyyyd, xxyyd, xxxxd = DSGE.var_approx_state_space(system[:TTT], system[:RRR], system[:QQ],
                                                 system[:DD], system[:ZZ], system[:EE],
                                                 zeros(size(system[:ZZ], 1),
                                                       DSGE.n_shocks_exogenous(m)),
                                                 4; get_covariances = true)
    yyyydc, xxyydc, xxxxdc = DSGE.var_approx_state_space(system[:TTT], system[:RRR], system[:QQ],
                                                         system[:DD], system[:ZZ], system[:EE],
                                                         zeros(size(system[:ZZ], 1),
                                                               DSGE.n_shocks_exogenous(m)),
                                                         4; get_covariances = true,
                                                         include_constant = true)
    β, Σ = DSGE.var_approx_state_space(system[:TTT], system[:RRR], system[:QQ], system[:DD],
                                  system[:ZZ], system[:EE],
                                  zeros(size(system[:ZZ], 1), DSGE.n_shocks_exogenous(m)),
                                  4; get_covariances = false)
    βc, Σc = DSGE.var_approx_state_space(system[:TTT], system[:RRR], system[:QQ], system[:DD],
                                  system[:ZZ], system[:EE],
                                  zeros(size(system[:ZZ], 1), DSGE.n_shocks_exogenous(m)),
                                  4; get_covariances = false, include_constant = true)

    expmat = matread("reference/exp_var_approx_state_space.mat")
    @test @test_matrix_approx_eq yyyyd expmat["yyyyd"]
    @test @test_matrix_approx_eq xxyyd expmat["xxyyd"][2:end, :]
    @test @test_matrix_approx_eq xxxxd expmat["xxxxd"][2:end, 2:end]
    @test @test_matrix_approx_eq yyyydc expmat["yyyyd"]
    @test @test_matrix_approx_eq xxyydc expmat["xxyyd"]
    @test @test_matrix_approx_eq xxxxdc expmat["xxxxd"]

    expβ = \(expmat["xxxxd"][2:end, 2:end], expmat["xxyyd"][2:end, :])
    expΣ = expmat["yyyyd"] - expmat["xxyyd"][2:end, :]' * expβ
    expΣ += expΣ'
    expΣ ./= 2
    @test @test_matrix_approx_eq β expβ
    @test @test_matrix_approx_eq Σ expΣ

    expβc = \(expmat["xxxxd"], expmat["xxyyd"])
    expΣc = expmat["yyyyd"] - expmat["xxyyd"]' * expβc
    expΣc += expΣc'
    expΣc ./= 2
    @test @test_matrix_approx_eq βc expβc
    @test @test_matrix_approx_eq Σc expΣc
end

nothing
