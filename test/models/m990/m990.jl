using DSGE
using HDF5, Base.Test, Distributions
include("../../util.jl")

path = dirname(@__FILE__)

### Model
model = Model990()

### Parameters

for θ in model.parameters
    !isa(θ,AbstractParameter) && error()
    θ.fixed && continue

    (left, right) = θ.valuebounds
    @test left < θ.value < right
end

### Model indices

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

### Equilibrium conditions
Γ0, Γ1, C, Ψ, Π = eqcond(model)

# Matrices are of expected dimensions
@test size(Γ0) == (60, 60)
@test size(Γ1) == (60, 60)
@test size(C) == (60, 1)
@test size(Ψ) == (60, 16)
@test size(Π) == (60, 13)

# Check output matrices against reference output (ϵ = 1e-4)
h5 = h5open("$path/eqcond.h5")
Γ0_ref = read(h5, "G0")
Γ1_ref = read(h5, "G1")
C_ref  = reshape(read(h5, "C"), 60, 1)
Ψ_ref  = read(h5, "PSI")
Π_ref  = read(h5, "PIE")
close(h5)

@test_matrix_approx_eq Γ0_ref Γ0
@test_matrix_approx_eq Γ1_ref Γ1
@test_matrix_approx_eq C_ref C
@test_matrix_approx_eq Ψ_ref Ψ
@test_matrix_approx_eq Π_ref Π

# ### Measurement equation
expect = Dict{Symbol, Matrix}()
h5 = h5open("$path/measurement.h5")
expect[:ZZ] = read(h5, "ZZ")
expect[:DD] = reshape(read(h5, "DD"), 12, 1)
expect[:QQ] = read(h5, "QQ")
expect[:EE] = read(h5, "EE")
expect[:MM]  = read(h5, "MM")
close(h5)

model = Model990()
TTT, RRR, CCC = solve(model)
actual = measurement(model, TTT, RRR, CCC)
for d in (:ZZ, :DD, :QQ, :EE, :MM)
    @test_matrix_approx_eq expect[d] actual[d]
end

nothing
