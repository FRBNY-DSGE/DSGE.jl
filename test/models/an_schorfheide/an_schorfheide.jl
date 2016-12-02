using DSGE
using HDF5, Base.Test, Distributions
include("../../util.jl")

path = dirname(@__FILE__)

### Model
model = AnSchorfheide()

### Model indices

# Endogenous states
endo = model.endogenous_states
@test length(endo) == 8
@test endo[:Eπ_t1] == 8

# Exogenous shocks
exo = model.exogenous_shocks
@test length(exo) == 3
@test exo[:R_sh] == 3

# Expectation shocks
ex = model.expected_shocks
@test length(ex) == 2
@test ex[:Eπ_sh] == 2

# Equations
eq = model.equilibrium_conditions
@test length(eq) == 8
@test eq[:shock_5] == 8

# Additional states
endo_new = model.endogenous_states_augmented
@test length(endo_new) == 0

# Observables
obs = model.observables
@test length(obs) == 3
@test obs[:obs_ffr] == 3

### Equilibrium conditions
Γ0, Γ1, C, Ψ, Π = eqcond(model)

# Matrices are of expected dimensions
@test size(Γ0) == (8, 8)
@test size(Γ1) == (8, 8)
@test size(C) == (8, 1)
@test size(Ψ) == (8, 3)
@test size(Π) == (8, 2)

nothing
