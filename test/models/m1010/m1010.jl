using DSGE
using Base.Test

path = dirname(@__FILE__)

### Model
model = Model1010("ss18")

### Parameters

for θ in model.parameters
    @test isa(θ, AbstractParameter)
    if !θ.fixed
        (left, right) = θ.valuebounds
        @test left < θ.value < right
    end
end


### Model indices

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

### Equilibrium conditions
Γ0, Γ1, C, Ψ, Π = eqcond(model)

# Transition and measurement equations
TTT, RRR, CCC = solve(model)
meas = measurement(model, TTT, RRR, CCC)

### Pseudo-measurement equation
pseudo, mapping = pseudo_measurement(model)

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
@test size(meas[:MM]) == (20,29)

@test size(TTT) == (91,91)
@test size(RRR) == (91,29)
@test size(CCC) == (91,)
