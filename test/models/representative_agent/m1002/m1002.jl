using DSGE
using Test

path = dirname(@__FILE__)

### Model
model = Model1002("ss2")

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
    @test length(endo) == 68
    @test endo[:Ez_t] == 60

    # Exogenous shocks
    exo = model.exogenous_shocks
    @test length(exo) == 24
    @test exo[:corepce_sh] == 16

    # Expectation shocks
    ex = model.expected_shocks
    @test length(ex) == 13
    @test ex[:EL_f_sh] == 12

    # Equations
    eq = model.equilibrium_conditions
    @test length(eq) == 68
    @test eq[:eq_Ez] == 59

    # Additional states
    endo_new = model.endogenous_states_augmented
    @test length(endo_new) == 16
    @test endo_new[:y_t1] == 69

    # Observables
    obs = model.observables
    @test length(obs) == 19
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
    @test size(Γ0) == (68, 68)
    @test size(Γ1) == (68, 68)
    @test size(C) == (68,)
    @test size(Ψ) == (68, 24)
    @test size(Π) == (68, 13)

    @test size(meas[:ZZ]) == (19,84)
    @test size(meas[:DD]) == (19,)
    @test size(meas[:QQ]) == (24,24)
    @test size(meas[:EE]) == (19,19)

    @test size(TTT) == (84,84)
    @test size(RRR) == (84,24)
    @test size(CCC) == (84,)
end
