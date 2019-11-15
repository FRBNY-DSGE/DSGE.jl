path = dirname(@__FILE__)

### Model
model = Model1010("ss20")

### Parameters

@testset "Test parameter bounds checking" begin
    for θ in model.parameters
        @test isa(θ, AbstractParameter)
        if !θ.fixed
            (left, right) = θ.valuebounds
            @test left < θ.value < right
        end
    end
end


### Model indices

@testset "Checking model indices" begin
    # Endogenous states
    endo = model.endogenous_states
    @test length(endo) == 73
    @test endo[:Ez_t] == 61

    # Exogenous shocks
    exo = model.exogenous_shocks
    @test length(exo) == 29
    @test exo[:corepce_sh] == 19

    # Expectation shocks
    ex = model.expected_shocks
    @test length(ex) == 13
    @test ex[:EL_f_sh] == 12

    # Equations
    eq = model.equilibrium_conditions
    @test length(eq) == 73
    @test eq[:eq_Ez] == 60

    # Additional states
    endo_new = model.endogenous_states_augmented
    @test length(endo_new) == 18
    @test endo_new[:y_t1] == 74

    # Observables
    obs = model.observables
    @test length(obs) == 20
    @test obs[:obs_tfp] == 12
end

### Equilibrium conditions
Γ0, Γ1, C, Ψ, Π = eqcond(model)

# Transition and measurement equations
TTT, RRR, CCC = solve(model)
meas = measurement(model, TTT, RRR, CCC)

### Pseudo-measurement equation
pseudo_meas = pseudo_measurement(model, TTT, RRR, CCC)

@testset "Check eqcond and state-space dimensions" begin
    # Matrices are of expected dimensions
    @test size(Γ0) == (73, 73)
    @test size(Γ1) == (73, 73)
    @test size(C) == (73,)
    @test size(Ψ) == (73, 29)
    @test size(Π) == (73, 13)

    @test size(meas[:ZZ]) == (20,91)
    @test size(meas[:DD]) == (20,)
    @test size(meas[:QQ]) == (29,29)
    @test size(meas[:EE]) == (20,20)

    @test size(TTT) == (91,91)
    @test size(RRR) == (91,29)
    @test size(CCC) == (91,)
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
