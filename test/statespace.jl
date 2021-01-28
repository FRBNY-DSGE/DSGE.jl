using DSGE, ModelConstructors, Dates, Test, LinearAlgebra, FileIO, Random, JLD2

writing_output = false # Write output for tests which use random values
if VERSION < v"1.5"
    ver = "111"
else
    ver = "150"
end

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

@testset "Using compute_system to update an existing system" begin

    # Check updating
    system1 = compute_system(m, system)
    system2 = compute_system(m, system; observables = collect(keys(m.observables))[1:end - 1],
                                  pseudo_observables = collect(keys(m.pseudo_observables))[1:end - 1],
                                  shocks = collect(keys(m.exogenous_shocks))[1:end - 1],
                                  states = collect(keys(m.endogenous_states))[1:end - 1])
    system3 = compute_system(m, system; zero_DD = true)
    system4 = compute_system(m, system; zero_DD_pseudo = true)

    @test system1[:TTT] ≈ system[:TTT]
    @test system1[:RRR] ≈ system[:RRR]
    @test system1[:CCC] ≈ system[:CCC]
    @test system1[:ZZ] ≈ system[:ZZ]
    @test system1[:DD] ≈ system[:DD]
    @test system1[:QQ] ≈ system[:QQ]
    @test system1[:EE] ≈ zeros(size(system[:ZZ], 1), size(system[:ZZ], 1))
    @test system1[:ZZ_pseudo] ≈ system[:ZZ_pseudo]
    @test system1[:DD_pseudo] ≈ system[:DD_pseudo]

    @test system2[:TTT] ≈ system[:TTT][1:end - 1, 1:end - 1]
    @test system2[:RRR] ≈ system[:RRR][1:end - 1, 1:end - 1]
    @test sum(abs.(system2[:CCC])) ≈ sum(abs.(system[:CCC][1:end - 1])) ≈ 0.
    @test system2[:ZZ] ≈ system[:ZZ][1:end - 1, 1:end - 1]
    @test system2[:DD] ≈ system[:DD][1:end - 1]
    @test system2[:QQ] ≈ system[:QQ][1:end - 1, 1:end - 1]
    @test system2[:EE] ≈ zeros(size(system[:ZZ], 1) - 1, size(system[:ZZ], 1) - 1)
    @test system2[:ZZ_pseudo] ≈ system[:ZZ_pseudo][1:end - 1, 1:end - 1]
    @test system2[:DD_pseudo] ≈ system[:DD_pseudo][1:end - 1]

    @test system3[:TTT] ≈ system[:TTT]
    @test system3[:RRR] ≈ system[:RRR]
    @test system3[:CCC] ≈ system[:CCC]
    @test system3[:ZZ] ≈ system[:ZZ]
    @test system3[:DD] ≈ zeros(length(system[:DD]))
    @test system3[:QQ] ≈ system[:QQ]
    @test system3[:EE] ≈ zeros(size(system[:ZZ], 1), size(system[:ZZ], 1))
    @test system3[:ZZ_pseudo] ≈ system[:ZZ_pseudo]
    @test system3[:DD_pseudo] ≈ system[:DD_pseudo]

    @test system4[:TTT] ≈ system[:TTT]
    @test system4[:RRR] ≈ system[:RRR]
    @test system4[:CCC] ≈ system[:CCC]
    @test system4[:ZZ] ≈ system[:ZZ]
    @test system4[:DD] ≈ system[:DD]
    @test system4[:QQ] ≈ system[:QQ]
    @test system4[:EE] ≈ zeros(size(system[:ZZ], 1), size(system[:ZZ], 1))
    @test system4[:ZZ_pseudo] ≈ system[:ZZ_pseudo]
    @test system4[:DD_pseudo] ≈ zeros(length(system[:DD_pseudo]))

    # Check errors
    @test_throws ErrorException compute_system(m, system; observables = [:blah])
    @test_throws ErrorException compute_system(m, system; states = [:blah])
    @test_throws KeyError compute_system(m, system; shocks = [:blah])
    @test_throws ErrorException compute_system(m, system; pseudo_observables = [:blah])
end

@testset "VAR approximation of state space" begin
    m = Model1002("ss10"; custom_settings =
                  Dict{Symbol,Setting}(:add_laborshare_measurement =>
                                       Setting(:add_laborshare_measurement, true),
                                       :add_NominalWageGrowth =>
                                       Setting(:add_NominalWageGrowth, true)))
    system = compute_system(m)
    system = compute_system(m, system; observables = [:obs_hours, :obs_gdpdeflator,
                                                      :laborshare_t, :NominalWageGrowth],
                            shocks = collect(keys(m.exogenous_shocks)))
    yyyyd, xxyyd, xxxxd = DSGE.var_approx_state_space(system[:TTT], system[:RRR], system[:QQ],
                                                 system[:DD], system[:ZZ], system[:EE],
                                                 zeros(size(system[:ZZ], 1),
                                                       DSGE.n_shocks_exogenous(m)),
                                                 4; get_population_moments = true)
    yyyydc, xxyydc, xxxxdc = DSGE.var_approx_state_space(system[:TTT], system[:RRR], system[:QQ],
                                                         system[:DD], system[:ZZ], system[:EE],
                                                         zeros(size(system[:ZZ], 1),
                                                               DSGE.n_shocks_exogenous(m)),
                                                         4; get_population_moments = true,
                                                         use_intercept = true)
    β, Σ = DSGE.var_approx_state_space(system[:TTT], system[:RRR], system[:QQ], system[:DD],
                                  system[:ZZ], system[:EE],
                                  zeros(size(system[:ZZ], 1), DSGE.n_shocks_exogenous(m)),
                                  4; get_population_moments = false)
    βc, Σc = DSGE.var_approx_state_space(system[:TTT], system[:RRR], system[:QQ], system[:DD],
                                  system[:ZZ], system[:EE],
                                  zeros(size(system[:ZZ], 1), DSGE.n_shocks_exogenous(m)),
                                  4; get_population_moments = false, use_intercept = true)

    expmat = load("reference/exp_var_approx_state_space.jld2")
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

    # Check DSGEVAR automates this properly
    dsgevar = DSGEVAR(m)
    DSGE.update!(dsgevar, shocks = collect(keys(m.exogenous_shocks)),
                 observables = [:obs_hours, :obs_gdpdeflator, :laborshare_t, :NominalWageGrowth],
                 lags = 4, λ = Inf)
    yyyyd, xxyyd, xxxxd = compute_system(dsgevar; get_population_moments = true)
    yyyydc, xxyydc, xxxxdc = compute_system(dsgevar; get_population_moments = true, use_intercept = true)
    β, Σ = compute_system(dsgevar)
    βc, Σc = compute_system(dsgevar; use_intercept = true)

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

@testset "VAR using DSGE as a prior" begin
    m = Model1002("ss10"; custom_settings =
                  Dict{Symbol,Setting}(:add_laborshare_measurement =>
                                       Setting(:add_laborshare_measurement, true),
                                       :add_NominalWageGrowth =>
                                       Setting(:add_NominalWageGrowth, true)))
    dsgevar = DSGEVAR(m)
    jlddata = load(joinpath(dirname(@__FILE__), "reference/test_dsgevar_lambda_irfs.jld2"))
    DSGE.update!(dsgevar, shocks = collect(keys(m.exogenous_shocks)),
                 observables = [:obs_hours, :obs_gdpdeflator, :laborshare_t, :NominalWageGrowth],
                 lags = 4, λ = Inf)
    data = jlddata["data"]
    yyyydc1, xxyydc1, xxxxdc1 = compute_system(dsgevar; get_population_moments = true, use_intercept = true)
    βc1, Σc1 = compute_system(dsgevar; use_intercept = true)

    yyyyd2, xxyyd2, xxxxd2 = compute_system(dsgevar, data; get_population_moments = true)
    β2, Σ2 = compute_system(dsgevar, data)

    # Check when λ = Inf
    @test @test_matrix_approx_eq yyyydc1 yyyyd2
    @test @test_matrix_approx_eq xxyydc1 xxyyd2
    @test @test_matrix_approx_eq xxxxdc1 xxxxd2
    @test @test_matrix_approx_eq βc1 β2
    @test @test_matrix_approx_eq Σc1 Σ2

    # Check when λ is finite
    DSGE.update!(dsgevar, λ = 1.)
    Random.seed!(1793) # need to seed for this
    yyyyd, xxyyd, xxxxd = compute_system(dsgevar, data; get_population_moments = true)
    β, Σ = compute_system(dsgevar, data)

    if writing_output
        JLD2.jldopen("reference/test_dsgevar_lambda_irfs_statespace_output_version=" * ver * ".jld2",
            true, true, true, IOStream) do file
            write(file, "exp_data_beta", β)
            write(file, "exp_data_sigma", Σ)
        end
    end

    file = jldopen("reference/test_dsgevar_lambda_irfs_statespace_output_version=" * ver * ".jld2", "r")
    saved_β = read(file, "exp_data_beta")
    saved_Σ = read(file, "exp_data_sigma")
    close(file)

    @test @test_matrix_approx_eq saved_β β
    @test @test_matrix_approx_eq saved_Σ Σ
end

@testset "VECM approximation of state space" begin
    matdata = load("reference/vecm_approx_state_space.jld2")
    nobs = Int(matdata["nvar"])
    p = Int(matdata["nlags"])
    coint = Int(matdata["coint"])
    TTT = matdata["TTT"]
    RRR = matdata["RRR"]
    ZZ = matdata["ZZ"]
    DD = vec(matdata["DD"])
    QQ = matdata["QQ"]
    EE = matdata["EE"]
    MM = matdata["MM"]
    yyyyd, xxyyd, xxxxd = DSGE.vecm_approx_state_space(TTT, RRR, QQ,
                                                       DD, ZZ, EE, MM, nobs,
                                                       p, coint; get_population_moments = true,
                                                       test_GA0 = matdata["GA0"])
    yyyydc, xxyydc, xxxxdc = DSGE.vecm_approx_state_space(TTT, RRR, QQ,
                                                          DD, ZZ, EE, MM, nobs,
                                                          p, coint; get_population_moments = true,
                                                          use_intercept = true, test_GA0 = matdata["GA0"])
    β, Σ = DSGE.vecm_approx_state_space(TTT, RRR, QQ, DD,
                                        ZZ, EE, MM, nobs, p, coint;
                                        get_population_moments = false,
                                        test_GA0 = matdata["GA0"])
    βc, Σc = DSGE.vecm_approx_state_space(TTT, RRR, QQ, DD,
                                          ZZ, EE, MM, nobs, p, coint;
                                          get_population_moments = false, use_intercept = true,
                                          test_GA0 = matdata["GA0"])

    no_int_inds = vcat(1:coint, coint + 2:size(matdata["xxyyd"], 1))
    @test @test_matrix_approx_eq yyyyd matdata["yyyyd"]
    @test @test_matrix_approx_eq xxyyd matdata["xxyyd"][no_int_inds, :]
    @test @test_matrix_approx_eq xxxxd matdata["xxxxd"][no_int_inds, no_int_inds]
    @test @test_matrix_approx_eq yyyydc matdata["yyyyd"]
    @test @test_matrix_approx_eq xxyydc matdata["xxyyd"]
    @test @test_matrix_approx_eq xxxxdc matdata["xxxxd"]

    expβ = \(matdata["xxxxd"][no_int_inds, no_int_inds], matdata["xxyyd"][no_int_inds, :])
    expΣ = matdata["yyyyd"] - matdata["xxyyd"][no_int_inds, :]' * expβ
    @test @test_matrix_approx_eq β expβ
    @test @test_matrix_approx_eq Σ expΣ

    expβc = matdata["xxxxd"] \ matdata["xxyyd"]
    expΣc = matdata["yyyyd"] - matdata["xxyyd"]' * expβc
    @test @test_matrix_approx_eq βc expβc
    @test @test_matrix_approx_eq Σc expΣc

    # Check DSGEVECM automates this properly. Need a proper model before being able to do this though.
    # m = Model1002("ss10")
    # dsgevecm = DSGEVECM(m)
    # DSGE.update!(dsgevecm, shocks = collect(keys(m.exogenous_shocks)),
    #              observables = [:obs_hours, :obs_gdpdeflator, :laborshare_t, :NominalWageGrowth],
    #              lags = 4, λ = Inf)
    # yyyyd, xxyyd, xxxxd = compute_system(dsgevecm; get_population_moments = true)
    # yyyydc, xxyydc, xxxxdc = compute_system(dsgevecm; get_population_moments = true, use_intercept = true)
    # β, Σ = compute_system(dsgevecm)
    # βc, Σc = compute_system(dsgevecm; use_intercept = true)

    # @test @test_matrix_approx_eq yyyyd matdata["yyyyd"]
    # @test @test_matrix_approx_eq xxyyd matdata["xxyyd"][2:end, :]
    # @test @test_matrix_approx_eq xxxxd matdata["xxxxd"][2:end, 2:end]
    # @test @test_matrix_approx_eq yyyydc matdata["yyyyd"]
    # @test @test_matrix_approx_eq xxyydc matdata["xxyyd"]
    # @test @test_matrix_approx_eq xxxxdc matdata["xxxxd"]

    # expβ = \(matdata["xxxxd"][2:end, 2:end], matdata["xxyyd"][2:end, :])
    # expΣ = matdata["yyyyd"] - matdata["xxyyd"][2:end, :]' * expβ
    # expΣ += expΣ'
    # expΣ ./= 2
    # @test @test_matrix_approx_eq β expβ
    # @test @test_matrix_approx_eq Σ expΣ

    # expβc = \(matdata["xxxxd"], matdata["xxyyd"])
    # expΣc = matdata["yyyyd"] - matdata["xxyyd"]' * expβc
    # expΣc += expΣc'
    # expΣc ./= 2
    # @test @test_matrix_approx_eq βc expβc
    # @test @test_matrix_approx_eq Σc expΣc
end

@testset "Updating a system for a DSGEVECM" begin
    dsge = AnSchorfheide()
    m = DSGEVECM(dsge)
    sys = compute_system(dsge)
    Dout1 = DSGE.compute_DD_coint_add(m, sys, [:obs_gdp, :obs_cpi])
    @info "The following 2 warnings about an empty vector are expected"
    sys_dsgevecm, Dout2 = compute_system(m, sys; get_DD_coint_add = true,
                                         cointegrating_add = [:obs_gdp, :obs_cpi])
    Dout3 = DSGE.compute_DD_coint_add(m, sys, Vector{Symbol}(undef, 0))
    _, Dout4 = compute_system(m, sys; get_DD_coint_add = true)

    # Check computing DD_coint_add
    @test Dout1 == sys[:DD][1:2]
    @test Dout2 == sys[:DD][1:2]
    @test isempty(Dout3)
    @test isempty(Dout4)

    # Check the system looks right
    @test @test_matrix_approx_eq sys[:TTT] sys_dsgevecm[:TTT]
    @test @test_matrix_approx_eq sys[:RRR] sys_dsgevecm[:RRR]
    @test @test_matrix_approx_eq sys[:CCC] sys_dsgevecm[:CCC]
    @test @test_matrix_approx_eq sys[:ZZ] sys_dsgevecm[:ZZ]
    @test @test_matrix_approx_eq sys[:DD] sys_dsgevecm[:DD]
    @test @test_matrix_approx_eq sys[:QQ] sys_dsgevecm[:QQ]
    @test @test_matrix_approx_eq sys_dsgevecm[:EE] zeros(size(sys_dsgevecm[:EE]))
end

@testset "RegimeSwitchingSystem" begin
    transitions         = Vector{Transition{Float64}}(undef, 3)
    measurements        = Vector{Measurement{Float64}}(undef, 3)
    pseudo_measurements = Vector{PseudoMeasurement{Float64}}(undef, 3)
    for i in 1:3
        transitions[i]         = system.transition
        measurements[i]        = system.measurement
        pseudo_measurements[i] = system.pseudo_measurement
    end

    regswitch_sys1 = RegimeSwitchingSystem(transitions, measurements, pseudo_measurements)
    regswitch_sys2 = RegimeSwitchingSystem(transitions, measurements)
    regswitch_sys3 = RegimeSwitchingSystem([system, system, system])
    regswitch_sys4 = copy(regswitch_sys1)
    regswitch_sys5 = deepcopy(regswitch_sys1)

    # Test access/utility functions first
    @test regswitch_sys1[:regimes]             == 1:3
    @test regswitch_sys1[:transitions]         == transitions
    @test regswitch_sys1[:measurements]        == measurements
    @test regswitch_sys1[:pseudo_measurements] == pseudo_measurements
    @test n_regimes(regswitch_sys1)            == 3
    for i in 1:3
        # Creating a System for a specific regime
        for tmpsys in [System(regswitch_sys1, i), System(regswitch_sys3, i),
                       System(regswitch_sys4, i), System(regswitch_sys5, i)]
            @test isa(tmpsys, System)
            @test @test_matrix_approx_eq tmpsys[:TTT] system[:TTT]
            @test @test_matrix_approx_eq tmpsys[:RRR] system[:RRR]
            @test @test_matrix_approx_eq tmpsys[:CCC] system[:CCC]
            @test @test_matrix_approx_eq tmpsys[:ZZ]  system[:ZZ]
            @test @test_matrix_approx_eq tmpsys[:DD]  system[:DD]
            @test @test_matrix_approx_eq tmpsys[:QQ]  system[:QQ]
            @test @test_matrix_approx_eq tmpsys[:EE]  system[:EE]
            @test @test_matrix_approx_eq tmpsys[:ZZ_pseudo]  system[:ZZ_pseudo]
            @test @test_matrix_approx_eq tmpsys[:DD_pseudo]  system[:DD_pseudo]
        end

        tmpsys2 = regswitch_sys2[i] # this creates a System underneath the hood
        @test isa(tmpsys2, System)
        @test @test_matrix_approx_eq tmpsys2[:TTT] system[:TTT]
        @test @test_matrix_approx_eq tmpsys2[:RRR] system[:RRR]
        @test @test_matrix_approx_eq tmpsys2[:CCC] system[:CCC]
        @test @test_matrix_approx_eq tmpsys2[:ZZ]  system[:ZZ]
        @test @test_matrix_approx_eq tmpsys2[:DD]  system[:DD]
        @test @test_matrix_approx_eq tmpsys2[:QQ]  system[:QQ]
        @test @test_matrix_approx_eq tmpsys2[:EE]  system[:EE]
        @test all(tmpsys2[:ZZ_pseudo] .== 0.)
        @test all(tmpsys2[:DD_pseudo] .== 0.)

        # Accessing specific data types or regime matrices
        @test @test_matrix_approx_eq regswitch_sys1[i, :transition][:TTT] transitions[i][:TTT]
        @test @test_matrix_approx_eq regswitch_sys1[i, :transition][:RRR] transitions[i][:RRR]
        @test @test_matrix_approx_eq regswitch_sys1[i, :transition][:CCC] transitions[i][:CCC]
        @test @test_matrix_approx_eq regswitch_sys1[i, :TTT]              transitions[i][:TTT]
        @test @test_matrix_approx_eq regswitch_sys1[i, :RRR]              transitions[i][:RRR]
        @test @test_matrix_approx_eq regswitch_sys1[i, :CCC]              transitions[i][:CCC]

        @test @test_matrix_approx_eq regswitch_sys1[i, :measurement][:ZZ] measurements[i][:ZZ]
        @test @test_matrix_approx_eq regswitch_sys1[i, :measurement][:DD] measurements[i][:DD]
        @test @test_matrix_approx_eq regswitch_sys1[i, :measurement][:QQ] measurements[i][:QQ]
        @test @test_matrix_approx_eq regswitch_sys1[i, :measurement][:EE] measurements[i][:EE]
        @test @test_matrix_approx_eq regswitch_sys1[i, :ZZ]               measurements[i][:ZZ]
        @test @test_matrix_approx_eq regswitch_sys1[i, :DD]               measurements[i][:DD]
        @test @test_matrix_approx_eq regswitch_sys1[i, :QQ]               measurements[i][:QQ]
        @test @test_matrix_approx_eq regswitch_sys1[i, :EE]               measurements[i][:EE]

        @test @test_matrix_approx_eq regswitch_sys1[i, :pseudo_measurement][:ZZ_pseudo] pseudo_measurements[i][:ZZ_pseudo]
        @test @test_matrix_approx_eq regswitch_sys1[i, :pseudo_measurement][:DD_pseudo] pseudo_measurements[i][:DD_pseudo]
        @test @test_matrix_approx_eq regswitch_sys1[i, :ZZ_pseudo]                      pseudo_measurements[i][:ZZ_pseudo]
        @test @test_matrix_approx_eq regswitch_sys1[i, :DD_pseudo]                      pseudo_measurements[i][:DD_pseudo]

        # Check errors when trying to access values
        @test_throws KeyError regswitch_sys1[1, :a]
        @test_throws KeyError regswitch_sys1[:a]
        @test_throws BoundsError regswitch_sys1[4, :transition]
        @test_throws BoundsError regswitch_sys1[4]
        @test_throws BoundsError System(regswitch_sys1, 4)

        # Check copying and deepcopying
        oldval = copy(regswitch_sys1[i, :TTT][1, 1])
        regswitch_sys1[i, :TTT][1, 1] = oldval + 1.
        @test @test_matrix_approx_eq regswitch_sys1[i, :TTT] regswitch_sys4[i, :TTT]
        @test !(regswitch_sys1[i, :TTT] ≈ regswitch_sys5[i, :TTT])
        regswitch_sys1[i, :TTT][1, 1] = oldval
    end
end

@testset "Implement alternative policy using regime_eqcond_info" begin
    output_vars = [:forecastobs, :histobs, :histpseudo, :forecastpseudo]

    m = Model1002("ss10", custom_settings = Dict{Symbol, Setting}(:add_altpolicy_pgap =>
                                                                  Setting(:add_altpolicy_pgap, true)))
    m <= Setting(:date_forecast_start, Date(2020, 6, 30))

    m <= Setting(:regime_switching, true)
    m <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m),
                                                  2 => Date(2020, 3, 31),
                                                  3 => Date(2020, 6, 30)))
    m = setup_regime_switching_inds!(m)
    m <= Setting(:forecast_horizons, 12)

    fp = dirname(@__FILE__)
    df = load(joinpath(fp, "reference", "regime_switch_data.jld2"), "regime_switch_df_none")

    m <= Setting(:replace_eqcond, true)
    m <= Setting(:regime_eqcond_info, Dict{Int, DSGE.EqcondEntry}(
                                                                  3 => DSGE.EqcondEntry(DSGE.ngdp(), [1., 0.])))
    m <= Setting(:pgap_type, :ngdp)
    m <= Setting(:pgap_value, 12.0)

    m <= Setting(:gensys2, false) # don't treat as a temporary policy
    sys1 = compute_system(m)
    #fcast_altperm = DSGE.forecast_one_draw(m, :mode, :full, output_vars, map(x -> x.value, m.parameters),
    #                                       df; regime_switching = true, n_regimes = get_setting(m, :n_regimes))

    # Testing ore basic permanent NGDP
    m = Model1002("ss10"; custom_settings = Dict{Symbol, Setting}(:add_altpolicy_pgap =>
                                                                    Setting(:add_altpolicy_pgap, true)))
    m <= Setting(:forecast_horizons, 12)
    m <= Setting(:date_forecast_start, Date(2020, 6, 30))
    m <= Setting(:regime_switching, true)
    m <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m),
                                                2 => Date(2020, 3, 31),
                                                3 => Date(2020, 6, 30)))
    m = setup_regime_switching_inds!(m)

    m <= Setting(:pgap_value, 12.0)
    m <= Setting(:pgap_type, :ngdp)
    m <= Setting(:regime_eqcond_info, Dict{Int, DSGE.EqcondEntry}(i => DSGE.EqcondEntry(DSGE.ngdp(), [1., 0.]) for i in 1:3))
    m <= Setting(:replace_eqcond, true)
    sys2 = compute_system(m)

    for i in 1:3
        @test !(sys1[1, :TTT] ≈ sys2[i, :TTT])
        @test !(sys1[2, :TTT] ≈ sys2[i, :TTT])
        @test sys1[3, :TTT] ≈ sys2[i, :TTT]
    end
    @test sys1[1, :TTT] ≈ sys1[2, :TTT]
end

@testset "Calculating transition matrices/vectors for k-periods ahead expectations and expected sums" begin
    # Test k-periods ahead expectations
    m = Model1002("ss10")
    sys_constant = compute_system(m)
    m <= Setting(:date_forecast_start, Date(2020, 9, 30))
    m <= Setting(:date_conditional_end, Date(2020, 9, 30))
    m <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m), 2 => Date(2020, 9, 30),
                                                3 => Date(2020, 12, 31), 4 => Date(2021, 3, 31),
                                                5 => Date(2021, 6, 30), 6 => Date(2021, 9, 30), 7 => Date(2021, 12, 31)))
    m <= Setting(:replace_eqcond, true)
    m <= Setting(:gensys2, true)
    m <= Setting(:regime_switching, true)
    regime_eqcond_info = Dict{Int, DSGE.EqcondEntry}()
    for i in 2:6
        regime_eqcond_info[i] = DSGE.EqcondEntry(DSGE.zero_rate(), [1., 0.])
    end
    regime_eqcond_info[7] = DSGE.EqcondEntry(AltPolicy(:historical, eqcond, solve), [1., 0.])
    m <= Setting(:regime_eqcond_info, regime_eqcond_info)
    setup_regime_switching_inds!(m)
    m <= Setting(:tvis_regime_eqcond_info, [regime_eqcond_info])
    m <= Setting(:tvis_information_set, [1:1, 2:7, 3:7, 4:7, 5:7, 6:7, 7:7])
    m <= Setting(:tvis_select_system, ones(Int, 7))
    sys = compute_system(m)

    T_acc1, C_acc1 = DSGE.k_periods_ahead_expectations(sys_constant[:TTT], sys_constant[:CCC], typeof(sys_constant[:TTT])[],
                                                       typeof(sys_constant[:CCC])[], 1, 3)
    T_acc2, C_acc2 = DSGE.k_periods_ahead_expectations(sys_constant[:TTT], sys_constant[:CCC], typeof(sys_constant[:TTT])[],
                                                       typeof(sys_constant[:CCC])[], 1, 3, 3)
    TTT_constant³ = sys_constant[:TTT] ^ 3
    @test T_acc1 ≈ TTT_constant³
    @test T_acc2 ≈ TTT_constant³
    @test C_acc1 ≈ sys_constant[:CCC]
    @test C_acc2 ≈ sys_constant[:CCC]

    T_acc3, C_acc3 = DSGE.k_periods_ahead_expectations(sys_constant[:TTT], sys_constant[:CCC] .+ .1, typeof(sys_constant[:TTT])[],
                                                       typeof(sys_constant[:CCC])[], 1, 3)
    TTT_constant³sum = (I - sys_constant[:TTT]) \ (I - TTT_constant³)
    @test T_acc3 ≈ TTT_constant³
    @test C_acc3 ≈ TTT_constant³sum * (sys_constant[:CCC] .+ .1)
    @test !(all(C_acc3 .≈ 0.))

    TTTs = [sys[i, :TTT] for i in 1:7]
    CCCs = [sys[i, :CCC] for i in 1:7]
    T_acc1, C_acc1 = DSGE.k_periods_ahead_expectations(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 2, 3)
    T_acc2, C_acc2 = DSGE.k_periods_ahead_expectations(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 2, 3, 7)
    T_acc3, C_acc3 = DSGE.k_periods_ahead_expectations(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 2, 5)
    T_acc4, C_acc4 = DSGE.k_periods_ahead_expectations(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 2, 7)
    T_acc5, C_acc5 = DSGE.k_periods_ahead_expectations(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 7, 3)

    @test T_acc1 ≈ sys[5, :TTT] * sys[4, :TTT] * sys[3, :TTT]
    @test T_acc2 ≈ sys[5, :TTT] * sys[4, :TTT] * sys[3, :TTT]
    @test T_acc3 ≈ sys[7, :TTT] * sys[6, :TTT] * sys[5, :TTT] * sys[4, :TTT] * sys[3, :TTT]
    @test T_acc4 ≈ sys[7, :TTT] * sys[7, :TTT] * sys[7, :TTT] * sys[6, :TTT] * sys[5, :TTT] * sys[4, :TTT] * sys[3, :TTT]
    @test T_acc5 ≈ sys[7, :TTT] ^ 3
    @test C_acc1 ≈ sys[5, :CCC] + sys[5, :TTT] * sys[4, :CCC] + sys[5, :TTT] * sys[4, :TTT] * sys[3, :CCC]
    @test C_acc2 ≈ sys[5, :CCC] + sys[5, :TTT] * sys[4, :CCC] + sys[5, :TTT] * sys[4, :TTT] * sys[3, :CCC]
    @test C_acc3 ≈ sys[7, :CCC] + sys[7, :TTT] * sys[6, :CCC] + sys[7, :TTT] * sys[6, :TTT] * sys[5, :CCC] +
        sys[7, :TTT] * sys[6, :TTT] * sys[5, :TTT] * sys[4, :CCC] + sys[7, :TTT] * sys[6, :TTT] * sys[5, :TTT] * sys[4, :TTT] * sys[3, :CCC]
    @test C_acc4 ≈ sys[7, :TTT]^2 * (sys[7, :CCC] + sys[7, :TTT] * sys[6, :CCC] + sys[7, :TTT] * sys[6, :TTT] * sys[5, :CCC] +
                                     sys[7, :TTT] * sys[6, :TTT] * sys[5, :TTT] * sys[4, :CCC] +
                                     sys[7, :TTT] * sys[6, :TTT] * sys[5, :TTT] * sys[4, :TTT] * sys[3, :CCC])
    @test all(C_acc5 .≈ 0.)

    CCCs = [sys[i, :CCC] .+ .1 for i in 1:7]
    T_acc1, C_acc1 = DSGE.k_periods_ahead_expectations(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 2, 3)
    T_acc2, C_acc2 = DSGE.k_periods_ahead_expectations(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 2, 3, 7)
    T_acc3, C_acc3 = DSGE.k_periods_ahead_expectations(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 2, 5)
    T_acc4, C_acc4 = DSGE.k_periods_ahead_expectations(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 2, 7)
    T_acc5, C_acc5 = DSGE.k_periods_ahead_expectations(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 7, 3)

    @test T_acc1 ≈ sys[5, :TTT] * sys[4, :TTT] * sys[3, :TTT]
    @test T_acc2 ≈ sys[5, :TTT] * sys[4, :TTT] * sys[3, :TTT]
    @test T_acc3 ≈ sys[7, :TTT] * sys[6, :TTT] * sys[5, :TTT] * sys[4, :TTT] * sys[3, :TTT]
    @test T_acc4 ≈ sys[7, :TTT] * sys[7, :TTT] * sys[7, :TTT] * sys[6, :TTT] * sys[5, :TTT] * sys[4, :TTT] * sys[3, :TTT]
    @test T_acc5 ≈ sys[7, :TTT] ^ 3
    @test C_acc1 ≈ CCCs[5] + sys[5, :TTT] * CCCs[4] + sys[5, :TTT] * sys[4, :TTT] * CCCs[3]
    @test C_acc2 ≈ CCCs[5] + sys[5, :TTT] * CCCs[4] + sys[5, :TTT] * sys[4, :TTT] * CCCs[3]
    @test C_acc3 ≈ CCCs[7] + sys[7, :TTT] * CCCs[6] + sys[7, :TTT] * sys[6, :TTT] * CCCs[5] +
        sys[7, :TTT] * sys[6, :TTT] * sys[5, :TTT] * CCCs[4] + sys[7, :TTT] * sys[6, :TTT] * sys[5, :TTT] * sys[4, :TTT] * CCCs[3]
    @test C_acc4 ≈ (I + sys[7, :TTT]) * CCCs[7] + sys[7, :TTT]^2 * (CCCs[7] + sys[7, :TTT] * CCCs[6] + sys[7, :TTT] * sys[6, :TTT] * CCCs[5] +
                                                                    sys[7, :TTT] * sys[6, :TTT] * sys[5, :TTT] * CCCs[4] +
                                                                    sys[7, :TTT] * sys[6, :TTT] * sys[5, :TTT] * sys[4, :TTT] * CCCs[3])
    @test !all(C_acc5 .≈ 0.)

    m = Model1002("ss10")
    m <= Setting(:tvis_information_set, [1:1, 2:7, 3:7, 4:7, 5:7, 6:7, 7:7])
    sys_constant = compute_system(m; tvis = true)
    m <= Setting(:date_forecast_start, Date(2020, 9, 30))
    m <= Setting(:date_conditional_end, Date(2020, 9, 30))
    m <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m), 2 => Date(2020, 9, 30),
                                                3 => Date(2020, 12, 31), 4 => Date(2021, 3, 31),
                                                5 => Date(2021, 6, 30), 6 => Date(2021, 9, 30), 7 => Date(2021, 12, 31)))
    m <= Setting(:replace_eqcond, true)
    m <= Setting(:gensys2, true)
    m <= Setting(:regime_switching, true)
    regime_eqcond_info = Dict{Int, DSGE.EqcondEntry}()
    for i in 2:6
        regime_eqcond_info[i] = DSGE.EqcondEntry(DSGE.zero_rate(), [1., 0.])
    end
    regime_eqcond_info[7] = DSGE.EqcondEntry(AltPolicy(:historical, eqcond, solve), [1., 0.])
    m <= Setting(:regime_eqcond_info, regime_eqcond_info)
    setup_regime_switching_inds!(m)
    sys = compute_system(m)

    T_acc1, C_acc1 = DSGE.k_periods_ahead_expected_sums(sys_constant[:TTT], sys_constant[:CCC], typeof(sys_constant[:TTT])[],
                                                        typeof(sys_constant[:CCC])[], 1, 3)
    T_acc2, C_acc2 = DSGE.k_periods_ahead_expected_sums(sys_constant[:TTT], sys_constant[:CCC], typeof(sys_constant[:TTT])[],
                                                        typeof(sys_constant[:CCC])[], 1, 3, 3)
    T_accum = Vector{typeof(sys_constant[:TTT])}(undef, 3)
    C_accum = Vector{typeof(sys_constant[:CCC])}(undef, 3)
    for i in 1:3
        T_accum[i], C_accum[i] = DSGE.k_periods_ahead_expectations(sys_constant[:TTT], sys_constant[:CCC], typeof(sys_constant[:TTT])[],
                                                                   typeof(sys_constant[:CCC])[], 1, i)
    end
    T_accum = sum(T_accum)
    C_accum = sum(C_accum)
    @test T_acc1 ≈ T_accum
    @test T_acc2 ≈ T_accum
    @test C_acc1 ≈ C_accum
    @test C_acc2 ≈ C_accum

    T_acc3, C_acc3 = DSGE.k_periods_ahead_expected_sums(sys_constant[:TTT], sys_constant[:CCC] .+ .1, typeof(sys_constant[:TTT])[],
                                                        typeof(sys_constant[:CCC])[], 1, 3)
    C_accum = Vector{typeof(sys_constant[:CCC])}(undef, 3)
    for i in 1:3
        _, C_accum[i] = DSGE.k_periods_ahead_expectations(sys_constant[:TTT], sys_constant[:CCC] .+ .1, typeof(sys_constant[:TTT])[],
                                                                   typeof(sys_constant[:CCC])[], 1, i)
    end
    @test T_acc3 ≈ T_accum
    @test C_acc3 ≈ sum(C_accum)

    TTTs = [sys[i, :TTT] for i in 1:7]
    CCCs = [sys[i, :CCC] for i in 1:7]
    T_acc1, C_acc1 = DSGE.k_periods_ahead_expected_sums(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 2, 3)
    T_acc2, C_acc2 = DSGE.k_periods_ahead_expected_sums(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 2, 3, 7)
    T_acc3, C_acc3 = DSGE.k_periods_ahead_expected_sums(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 2, 5)
    T_acc4, C_acc4 = DSGE.k_periods_ahead_expected_sums(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 2, 7)
    T_acc5, C_acc5 = DSGE.k_periods_ahead_expected_sums(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 7, 3)

    T_accum = Vector{typeof(sys[1, :TTT])}(undef, 3)
    C_accum = Vector{typeof(sys[1, :CCC])}(undef, 3)
    for i in 1:3
        T_accum[i], C_accum[i] = DSGE.k_periods_ahead_expectations(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 2, i)
    end
    @test T_acc1 ≈ sum(T_accum)
    @test T_acc2 ≈ sum(T_accum)
    @test C_acc1 ≈ sum(C_accum)
    @test C_acc2 ≈ sum(C_accum)

    T_accum = Vector{typeof(sys[1, :TTT])}(undef, 5)
    C_accum = Vector{typeof(sys[1, :CCC])}(undef, 5)
    for i in 1:5
        T_accum[i], C_accum[i] = DSGE.k_periods_ahead_expectations(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 2, i)
    end
    @test T_acc3 ≈ sum(T_accum)
    @test C_acc3 ≈ sum(C_accum)

    T_accum = Vector{typeof(sys[1, :TTT])}(undef, 7)
    C_accum = Vector{typeof(sys[1, :CCC])}(undef, 7)
    for i in 1:7
        T_accum[i], C_accum[i] = DSGE.k_periods_ahead_expectations(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 2, i)
    end
    @test T_acc4 ≈ sum(T_accum)
    @test C_acc4 ≈ sum(C_accum)

    T_accum = Vector{typeof(sys[1, :TTT])}(undef, 3)
    C_accum = Vector{typeof(sys[1, :CCC])}(undef, 3)
    for i in 1:3
        T_accum[i], C_accum[i] = DSGE.k_periods_ahead_expectations(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 7, i)
    end
    @test T_acc5 ≈ sum(T_accum)
    @test C_acc5 ≈ sum(C_accum)

    # Test k-periods ahead expected sum
    m = Model1002("ss10")
    sys_constant = compute_system(m)
    m <= Setting(:date_forecast_start, Date(2020, 9, 30))
    m <= Setting(:date_conditional_end, Date(2020, 9, 30))
    m <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m), 2 => Date(2020, 9, 30),
                                                3 => Date(2020, 12, 31), 4 => Date(2021, 3, 31),
                                                5 => Date(2021, 6, 30), 6 => Date(2021, 9, 30), 7 => Date(2021, 12, 31)))
    m <= Setting(:replace_eqcond, true)
    m <= Setting(:gensys2, true)
    m <= Setting(:regime_switching, true)
    regime_eqcond_info = Dict()
    for i in 2:6
        regime_eqcond_info[i] = DSGE.EqcondEntry(DSGE.zero_rate(), [1., 0.])
    end
    regime_eqcond_info[7] = DSGE.EqcondEntry(AltPolicy(:historical, eqcond, solve), [1., 0.])
    m <= Setting(:regime_eqcond_info, regime_eqcond_info)
    setup_regime_switching_inds!(m)
    sys = compute_system(m)

    T_acc1, C_acc1 = DSGE.k_periods_ahead_expected_sums(sys_constant[:TTT], sys_constant[:CCC], typeof(sys_constant[:TTT])[],
                                                        typeof(sys_constant[:CCC])[], 1, 3)
    T_acc2, C_acc2 = DSGE.k_periods_ahead_expected_sums(sys_constant[:TTT], sys_constant[:CCC], typeof(sys_constant[:TTT])[],
                                                        typeof(sys_constant[:CCC])[], 1, 3, 3)
    T_accum = Vector{typeof(sys_constant[:TTT])}(undef, 3)
    C_accum = Vector{typeof(sys_constant[:CCC])}(undef, 3)
    for i in 1:3
        T_accum[i], C_accum[i] = DSGE.k_periods_ahead_expectations(sys_constant[:TTT], sys_constant[:CCC], typeof(sys_constant[:TTT])[],
                                                                   typeof(sys_constant[:CCC])[], 1, i)
    end
    T_accum = sum(T_accum)
    C_accum = sum(C_accum)
    @test T_acc1 ≈ T_accum
    @test T_acc2 ≈ T_accum
    @test C_acc1 ≈ C_accum
    @test C_acc2 ≈ C_accum

    T_acc3, C_acc3 = DSGE.k_periods_ahead_expected_sums(sys_constant[:TTT], sys_constant[:CCC] .+ .1, typeof(sys_constant[:TTT])[],
                                                        typeof(sys_constant[:CCC])[], 1, 3)
    C_accum = Vector{typeof(sys_constant[:CCC])}(undef, 3)
    for i in 1:3
        _, C_accum[i] = DSGE.k_periods_ahead_expectations(sys_constant[:TTT], sys_constant[:CCC] .+ .1, typeof(sys_constant[:TTT])[],
                                                                   typeof(sys_constant[:CCC])[], 1, i)
    end
    @test T_acc3 ≈ T_accum
    @test C_acc3 ≈ sum(C_accum)

    TTTs = [sys[i, :TTT] for i in 1:7]
    CCCs = [sys[i, :CCC] for i in 1:7]
    T_acc1, C_acc1 = DSGE.k_periods_ahead_expected_sums(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 2, 3)
    T_acc2, C_acc2 = DSGE.k_periods_ahead_expected_sums(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 2, 3, 7)
    T_acc3, C_acc3 = DSGE.k_periods_ahead_expected_sums(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 2, 5)
    T_acc4, C_acc4 = DSGE.k_periods_ahead_expected_sums(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 2, 7)
    T_acc5, C_acc5 = DSGE.k_periods_ahead_expected_sums(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 7, 3)

    T_accum = Vector{typeof(sys[1, :TTT])}(undef, 3)
    C_accum = Vector{typeof(sys[1, :CCC])}(undef, 3)
    for i in 1:3
        T_accum[i], C_accum[i] = DSGE.k_periods_ahead_expectations(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 2, i)
    end
    @test T_acc1 ≈ sum(T_accum)
    @test T_acc2 ≈ sum(T_accum)
    @test C_acc1 ≈ sum(C_accum)
    @test C_acc2 ≈ sum(C_accum)

    T_accum = Vector{typeof(sys[1, :TTT])}(undef, 5)
    C_accum = Vector{typeof(sys[1, :CCC])}(undef, 5)
    for i in 1:5
        T_accum[i], C_accum[i] = DSGE.k_periods_ahead_expectations(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 2, i)
    end
    @test T_acc3 ≈ sum(T_accum)
    @test C_acc3 ≈ sum(C_accum)

    T_accum = Vector{typeof(sys[1, :TTT])}(undef, 7)
    C_accum = Vector{typeof(sys[1, :CCC])}(undef, 7)
    for i in 1:7
        T_accum[i], C_accum[i] = DSGE.k_periods_ahead_expectations(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 2, i)
    end
    @test T_acc4 ≈ sum(T_accum)
    @test C_acc4 ≈ sum(C_accum)

    T_accum = Vector{typeof(sys[1, :TTT])}(undef, 3)
    C_accum = Vector{typeof(sys[1, :CCC])}(undef, 3)
    for i in 1:3
        T_accum[i], C_accum[i] = DSGE.k_periods_ahead_expectations(sys[1, :TTT], sys[1, :CCC], TTTs, CCCs, 7, i)
    end
    @test T_acc5 ≈ sum(T_accum)
    @test C_acc5 ≈ sum(C_accum)
end

@testset "Time-Varying Information Set state space system" begin
    m = Model1002("ss10")
    system = compute_system(m)
    m <= Setting(:date_forecast_start, Date(2020, 9, 30))
    m <= Setting(:date_conditional_end, Date(2020, 9, 30))
    m <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m), 2 => Date(2020, 9, 30),
                                                3 => Date(2020, 12, 31), 4 => Date(2021, 3, 31),
                                                5 => Date(2021, 6, 30), 6 => Date(2021, 9, 30), 7 => Date(2021, 12, 31)))
    setup_regime_switching_inds!(m)
    m <= Setting(:replace_eqcond, true)
    m <= Setting(:gensys2, true)
    m <= Setting(:regime_switching, true)

    # Check that using the same regime_eqcond_info but different tvis_select_system yields the same results
    regime_eqcond_info1 = Dict()
    for i in 2:6
        regime_eqcond_info1[i] = DSGE.EqcondEntry(DSGE.zero_rate())#, [1., 0.])
    end
    regime_eqcond_info1[7] = DSGE.EqcondEntry(AltPolicy(:historical, eqcond, solve))#, [1., 0.])
    regime_eqcond_info2 = deepcopy(regime_eqcond_info1)
    m <= Setting(:tvis_regime_eqcond_info, [regime_eqcond_info1, regime_eqcond_info2])
    m <= Setting(:tvis_information_set, [1:1, 2:7, 3:7, 4:7, 5:7, 6:7, 7:7])
    m <= Setting(:tvis_select_system, ones(Int, 7))
    sys1 = DSGE.compute_tvis_system(m)
    m <= Setting(:tvis_select_system, fill(2, 7))
    sys2 = DSGE.compute_tvis_system(m)

    @test isa(sys1[:transitions], Vector{Vector{Transition{Float64}}})
    @test isa(sys1[:measurements], Vector{Measurement{Float64}})
    @test isa(sys1[:pseudo_measurements], Vector{PseudoMeasurement{Float64}})
    @test sys1[:information_set] == sys2.information_set
    @test all(sys1[:select] .!= sys2.select)
    @test sys1[1, :measurement][:ZZ] == sys1.measurements[1][:ZZ]
    @test sys1[1, :pseudo_measurement][:ZZ_pseudo] == sys1.pseudo_measurements[1][:ZZ_pseudo]
    @test sys1[1, :ZZ] == sys1.measurements[1][:ZZ]
    @test sys1[1, :DD] == sys1.measurements[1][:DD]
    @test sys1[1, :ZZ_pseudo] == sys1.pseudo_measurements[1][:ZZ_pseudo]
    @test sys1[1, 2][:TTT] == sys1.transitions[1][2][:TTT]
    @test sys1[1, 2, :TTT] == sys1.transitions[1][2][:TTT]
    @test n_regimes(sys1) == 7
    @test sys1[:regimes] == 1:7
    for reg in 1:n_regimes(sys1)
        sys1[reg, :ZZ] == sys2[reg, :ZZ]
    end

    # Check using the different regime_eqcond_info
    regime_eqcond_info2[6] = DSGE.EqcondEntry(AltPolicy(:historical, eqcond, solve))#, [1., 0.]) # shorten the ZLB
    m <= Setting(:tvis_regime_eqcond_info, [regime_eqcond_info1, regime_eqcond_info2])
    m <= Setting(:tvis_information_set, [1:1, 2:7, 3:7, 4:7, 5:7, 6:7, 7:7])
    m <= Setting(:tvis_select_system, [1, 1, 1, 1, 2, 2, 2])
    sys1 = DSGE.compute_tvis_system(m) # Does not seem to be doing regime-switching
    m <= Setting(:regime_eqcond_info, regime_eqcond_info1)
    m <= Setting(:tvis_regime_eqcond_info, [regime_eqcond_info1])
    m <= Setting(:tvis_select_system, ones(Int, 7))
    sys2 = compute_system(m; tvis = true)
    m <= Setting(:regime_eqcond_info, regime_eqcond_info2)
    m <= Setting(:tvis_regime_eqcond_info, [regime_eqcond_info2])
    m <= Setting(:tvis_select_system, ones(Int, 7)) # ones here b/c there's only one system for the TVIS to use
    sys3 = compute_system(m; tvis = true)
    for reg in 1:4
        @test sys1[reg, :ZZ] ≈ sys2[reg, :ZZ]
    end
    for reg in 5:7
        @test sys1[reg, :ZZ] ≈ sys3[reg, :ZZ]
    end
    for reg in 2:5
        @test sys2[reg, :ZZ] != sys3[reg, :ZZ]
    end
    for reg in 6:7 # same for 6 and 7 b/c in both cases, next period expectation does not use zero rate rule
        @test sys2[reg, :ZZ] == sys3[reg, :ZZ]
    end

    m <= Setting(:tvis_regime_eqcond_info, [regime_eqcond_info1, regime_eqcond_info2])
    m <= Setting(:tvis_select_system, [1, 1, 1, 1, 2, 2, 2])
    sys4 = compute_system(m; tvis = true)
    for reg in 1:4
        @test sys4[reg, :ZZ] == sys1[reg, :ZZ]
        @test sys4[reg, :ZZ] == sys2[reg, :ZZ]
        @test sys4[reg, :TTT] == sys1[1, reg, :TTT]
        @test sys4[reg, :TTT] == sys2[reg, :TTT]
    end
    for reg in 5:7
        @test sys4[reg, :ZZ] == sys1[reg, :ZZ]
        @test sys4[reg, :ZZ] == sys3[reg, :ZZ]
        @test sys4[reg, :TTT] == sys1[2, reg, :TTT]
        @test sys4[reg, :TTT] == sys3[reg, :TTT]
    end
    for reg in 2:4
        @test sys4[reg, :ZZ] != sys3[reg, :ZZ]
    end
    for reg in 5:5 # only to 5 b/c the difference in info sets about ZLB length
        @test sys4[reg, :ZZ] != sys2[reg, :ZZ]
    end
end

@testset "Using gensys2 in historical regimes" begin
    function set_regime_vals_fnct!(m, n)
        if n > 4
            for p in m.parameters
                if !isempty(p.regimes) && haskey(p.regimes, :value)
                    for i in 5:n
                        ModelConstructors.set_regime_val!(p, i, ModelConstructors.regime_val(p, 4))
                    end
                end
            end
        end
    end

    imperfect_cred_new = 1.
    imperfect_cred_old = 1. - imperfect_cred_new
    custom_set = Dict{Symbol,Setting}(:n_mon_anticipated_shocks =>
                                      Setting(:n_mon_anticipated_shocks, 6, "Number of anticipated policy shocks"),
                                      :imperfect_awareness_weights => Setting(:imperfect_awareness_weights,
                                                                              [imperfect_cred_new, imperfect_cred_old]),
                                      :alternative_policies => Setting(:alternative_policies, [DSGE.taylor_rule()]),
                                      :flexible_ait_policy_change => Setting(:flexible_ait_policy_change, true))

    m = Model1002("ss59", custom_settings = custom_set)
    usual_model_settings!(m, "201117", fcast_date = Date(2020, 9, 30))

    set_regime_vals_fnct!(m, 4 + 4)
    m <= Setting(:replace_eqcond, true)
    m <= Setting(:gensys2, true)
    reg_dates = deepcopy(get_setting(m, :regime_dates))
    regime_eqcond_info = Dict{Int, DSGE.EqcondEntry}()
    for (regind, date) in zip(4:(4 + 4), # 4 + 1 b/c zero for 4 periods and liftoff, add 3 to make correct regime
                              DSGE.quarter_range(reg_dates[4], DSGE.iterate_quarters(reg_dates[4], 4)))
        reg_dates[regind] = date
        if regind != 4 + 4
            regime_eqcond_info[regind] = DSGE.EqcondEntry(DSGE.zero_rate(), [1., 0.])
        end
    end
    regime_eqcond_info[8] = DSGE.EqcondEntry(DSGE.flexible_ait(), [imperfect_cred_new, imperfect_cred_old])
    m <= Setting(:regime_dates, reg_dates)
    m <= Setting(:regime_eqcond_info, regime_eqcond_info)
    setup_regime_switching_inds!(m; cond_type = :full, temp_altpolicy_in_cond_regimes = true)
    m <= Setting(:tvis_information_set, [1:1, 2:2, 3:3, [i:get_setting(m, :n_regimes) for i in 4:get_setting(m, :n_regimes)]...])

    sys_true = compute_system(m; tvis = true)
    m <= Setting(:uncertain_temp_altpol, true)
    sys_true_unc = compute_system(m; tvis = true)
    m <= Setting(:uncertain_temp_altpol, false)

    m <= Setting(:date_forecast_start, Date(2020, 12, 31))
    m <= Setting(:date_conditional_end, Date(2020, 12, 31))
    setup_regime_switching_inds!(m; cond_type = :full, temp_altpolicy_in_cond_regimes = true)
    sys = compute_system(m; tvis = true)
    for i in 1:get_setting(m, :n_regimes)
        @test sys_true[i, :TTT] ≈ sys[i, :TTT]
        @test sys_true[i, :ZZ] ≈ sys[i, :ZZ]
    end

    m <= Setting(:uncertain_temp_altpol, true)
    sys = compute_system(m; tvis = true)
    for i in 1:get_setting(m, :n_regimes)
        @test sys_true_unc[i, :TTT] ≈ sys[i, :TTT]
        @test sys_true_unc[i, :ZZ] ≈ sys[i, :ZZ]
    end
end

@testset "Measurement equation of forward-looking variables" begin
    function set_regime_vals_fnct!(m, n)
        if n > 4
            for p in m.parameters
                if !isempty(p.regimes) && haskey(p.regimes, :value)
                    for i in 5:n
                        ModelConstructors.set_regime_val!(p, i, ModelConstructors.regime_val(p, 4))
                    end
                end
            end
        end
    end
    imperfect_cred_new = .5
    imperfect_cred_old = 1. - imperfect_cred_new
    custom_set = Dict{Symbol,Setting}(:n_mon_anticipated_shocks =>
                                      Setting(:n_mon_anticipated_shocks, 6, "Number of anticipated policy shocks"),
                                      :imperfect_awareness_weights => Setting(:imperfect_awareness_weights,
                                                                              [imperfect_cred_new, imperfect_cred_old]),
                                      :alternative_policies => Setting(:alternative_policies, [DSGE.taylor_rule()]),
                                      :flexible_ait_policy_change => Setting(:flexible_ait_policy_change, true))

    m = Model1002("ss59", custom_settings = custom_set)
    usual_model_settings!(m, "201117", fcast_date = Date(2020, 9, 30))
    m <= Setting(:flexible_ait_policy_change, false) # was true just so initialize a bunch of settings

    set_regime_vals_fnct!(m, 4 + 4)
    m <= Setting(:replace_eqcond, true)
    m <= Setting(:gensys2, true)
    reg_dates = deepcopy(get_setting(m, :regime_dates))
    regime_eqcond_info = Dict{Int, DSGE.EqcondEntry}()
    for (regind, date) in zip(4:(4 + 4), # 4 + 1 b/c zero for 4 periods and liftoff, add 3 to make correct regime
                              DSGE.quarter_range(reg_dates[4], DSGE.iterate_quarters(reg_dates[4], 4)))
        reg_dates[regind] = date
        if regind != 4 + 4
            regime_eqcond_info[regind] = DSGE.EqcondEntry(DSGE.zero_rate(), [imperfect_cred_new, imperfect_cred_old])
        end
    end
    regime_eqcond_info[8] = DSGE.EqcondEntry(DSGE.flexible_ait(), [imperfect_cred_new, imperfect_cred_old])
    regime_eqcond_info_cp = deepcopy(regime_eqcond_info)
    m <= Setting(:regime_dates, reg_dates)
    m <= Setting(:regime_eqcond_info, regime_eqcond_info)
    setup_regime_switching_inds!(m; cond_type = :full, temp_altpolicy_in_cond_regimes = true)
# GET RID OF THIS I THINK    m <= Setting(:tvis_information_set, [1:1, 2:2, 3:3, [i:get_setting(m, :n_regimes) for i in 4:get_setting(m, :n_regimes)]...])

    m <= Setting(:regime_switching, false)
    m <= Setting(:uncertain_temp_altpol, false)
    m <= Setting(:uncertain_altpolicy, false)
    delete!(m.settings, :regime_eqcond_info)
    m <= Setting(:gensys2, false)
    sys_taylor = compute_system(m)

    m <= Setting(:regime_switching, true)
    m <= Setting(:uncertain_temp_altpol, false)
    m <= Setting(:regime_eqcond_info, deepcopy(regime_eqcond_info_cp))
    m <= Setting(:uncertain_altpolicy, false)
    m <= Setting(:gensys2, true)
    sys_perfcred = compute_system(m)

    m <= Setting(:regime_eqcond_info, deepcopy(regime_eqcond_info_cp))
    m <= Setting(:uncertain_temp_altpol, true)
    m <= Setting(:uncertain_altpolicy, true)
    sys_unczlb_uncalt = compute_system(m)

    for i in sort!(collect(keys(get_setting(m, :regime_eqcond_info))))
        @test sys_unczlb_uncalt[i, :ZZ] ≈ imperfect_cred_new * sys_perfcred[i, :ZZ] +
            imperfect_cred_old * sys_taylor[:ZZ]
        @test sys_unczlb_uncalt[i, :ZZ_pseudo] ≈ imperfect_cred_new * sys_perfcred[i, :ZZ_pseudo] +
            imperfect_cred_old * sys_taylor[:ZZ_pseudo]
        @test sys_unczlb_uncalt[i, :DD] ≈ imperfect_cred_new * sys_perfcred[i, :DD] +
            imperfect_cred_old * sys_taylor[:DD]
        @test sys_unczlb_uncalt[i, :DD_pseudo] ≈ imperfect_cred_new * sys_perfcred[i, :DD_pseudo] +
            imperfect_cred_old * sys_taylor[:DD_pseudo]
    end
end

nothing
