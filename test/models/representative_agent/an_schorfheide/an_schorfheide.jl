using DSGE
using Test

path = dirname(@__FILE__)

### Model
model = AnSchorfheide()

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
    @test length(endo) == 8
    @test endo[:z_t] == 6

    # Exogenous shocks
    exo = model.exogenous_shocks
    @test length(exo) == 3
    @test exo[:g_sh] == 2

    # Expectation shocks
    ex = model.expected_shocks
    @test length(ex) == 2
    @test ex[:Ey_sh] == 1

    # Equations
    eq = model.equilibrium_conditions
    @test length(eq) == 8
    @test eq[:eq_mp] == 3

    # Additional states
    endo_new = model.endogenous_states_augmented
    @test length(endo_new) == 0

    # Observables
    obs = model.observables
    @test length(obs) == 3
    @test obs[:obs_nominalrate] == 3
end

### Equilibrium conditions
Γ0, Γ1, C, Ψ, Π = eqcond(model)

# Transition and measurement equations
TTT, RRR, CCC = solve(model)
meas = measurement(model, TTT, RRR, CCC)

@testset "Check eqcond and state-space dimensions" begin
    # Matrices are of expected dimensions
    @test size(Γ0) == (8, 8)
    @test size(Γ1) == (8, 8)
    @test size(C) == (8,)
    @test size(Ψ) == (8, 3)
    @test size(Π) == (8, 2)

    @test size(meas[:ZZ]) == (3, 8)
    @test size(meas[:DD]) == (3,)
    @test size(meas[:QQ]) == (3, 3)
    @test size(meas[:EE]) == (3, 3)

    @test size(TTT) == (8,8)
    @test size(RRR) == (8,3)
    @test size(CCC) == (8,)
end
