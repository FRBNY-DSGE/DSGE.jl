path = dirname(@__FILE__)

### Model
sub = "ss3"
model = Model990(sub)

### Parameters

@testset "Test parameter bounds checking" begin
    for θ in model.parameters
        !isa(θ,AbstractParameter) && error()
        θ.fixed && continue

        (left, right) = θ.valuebounds
        @test left < θ.value < right
    end
end

### Model indices

@testset "Checking model indices" begin
    # Endogenous states
    endo = model.endogenous_states
    @test length(endo) == 60
    @test endo[:Ez_t] == 60

    # Exogenous shocks
    exo = model.exogenous_shocks
    @test length(exo) == 16
    @test exo[:corepce_sh] == 16

    # Expectation shocks
    ex = model.expected_shocks
    @test length(ex) == 13
    @test ex[:Erk_f_sh] == 13

    # Equations
    eq = model.equilibrium_conditions
    @test length(eq) == 60
    @test eq[:eq_Ez] == 60

    # Additional states
    endo_new = model.endogenous_states_augmented
    @test length(endo_new) == 12
    @test endo_new[:y_t1] == 61

    # Observables
    obs = model.observables
    @test length(obs) == 12
    @test obs[:obs_tfp] == 12
end

### Equilibrium conditions
Γ0, Γ1, C, Ψ, Π = eqcond(model)

# Matrices are of expected dimensions
@testset "Check eqcond and state-space dimensions" begin
    @test size(Γ0) == (60, 60)
    @test size(Γ1) == (60, 60)
    @test size(C) == (60,)
    @test size(Ψ) == (60, 16)
    @test size(Π) == (60, 13)
end

### Custom settings
custom_settings = Dict{Symbol, Setting}(
    :n_mon_anticipated_shocks => Setting(:n_mon_anticipated_shocks, 6))
model = Model990(sub, custom_settings = custom_settings)
@testset "Checking n_ant_shocks == 6" begin
    @test get_setting(model, :n_mon_anticipated_shocks) == 6
end

# Indices initialized correctly under custom settings

@testset "Checking model indices part 2" begin
    # Endogenous states
    endo = model.endogenous_states
    @test length(endo) == 66
    @test endo[:rm_tl6] == 66

    # Exogenous shocks
    exo = model.exogenous_shocks
    @test length(exo) == 22
    @test exo[:rm_shl6] == 22

    # Expectation shocks
    ex = model.expected_shocks
    @test length(ex) == 13
    @test ex[:Erk_f_sh] == 13

    # Equations
    eq = model.equilibrium_conditions
    @test length(eq) == 66
    @test eq[:eq_rml6] == 66

    # Additional states
    endo_new = model.endogenous_states_augmented
    @test length(endo_new) == 12
    @test endo_new[:y_t1] == 67

    # Observables
    obs = model.observables
    @test length(obs) == 18
    @test obs[:obs_nominalrate6] == 18
end

### Equilibrium conditions
Γ0, Γ1, C, Ψ, Π = eqcond(model)

@testset "Check eqcond and state-space dimensions part 2" begin
    # Matrices are of expected dimensions
    @test size(Γ0) == (66, 66)
    @test size(Γ1) == (66, 66)
    @test size(C) == (66,)
    @test size(Ψ) == (66, 22)
    @test size(Π) == (66, 13)
end

### Reference matrices
G0_exp, G1_exp, C_exp, PSI_exp, PIE_exp = h5open("$path/eqcond.h5") do file
    read(file, "G0"), read(file, "G1"), read(file, "C"),
    read(file, "PSI"), read(file, "PIE")
end
TTT_exp, RRR_exp, CCC_exp = h5open("$path/transition.h5") do file
    read(file, "TTT"), read(file, "RRR"), read(file, "CCC")
end
Z_exp, D_exp, Q_exp, E_exp = h5open("$path/measurement.h5") do file
    read(file, "ZZ"), read(file, "DD"), read(file, "QQ"), read(file, "EE")
end
Z_pseudo_exp, D_pseudo_exp = h5open("$path/pseudo_measurement.h5") do file
    read(file, "ZZ_pseudo"), read(file, "DD_pseudo")
end

Γ0, Γ1, C, Ψ, Π = eqcond(model)
TTT, RRR, CCC = solve(model)
meas = measurement(model, TTT, RRR, CCC)
pseudo_meas = pseudo_measurement(model, TTT, RRR, CCC)

@testset "Compare eqcond, transition, measurement, and pseudo-measurement matrices against reference matrices" begin
    @test @test_matrix_approx_eq Γ0 G0_exp
    @test @test_matrix_approx_eq Γ1 G1_exp
    @test @test_matrix_approx_eq C C_exp
    @test @test_matrix_approx_eq Ψ PSI_exp
    @test @test_matrix_approx_eq Π PIE_exp

    @test @test_matrix_approx_eq TTT TTT_exp
    @test @test_matrix_approx_eq RRR RRR_exp
    @test @test_matrix_approx_eq CCC CCC_exp

    @test @test_matrix_approx_eq meas[:ZZ] Z_exp
    @test @test_matrix_approx_eq meas[:DD] D_exp
    @test @test_matrix_approx_eq meas[:QQ] Q_exp
    @test @test_matrix_approx_eq meas[:EE] E_exp

    @test @test_matrix_approx_eq pseudo_meas[:ZZ_pseudo] Z_pseudo_exp
    @test @test_matrix_approx_eq pseudo_meas[:DD_pseudo] D_pseudo_exp
end


nothing
