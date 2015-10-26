using Distributions,  Base.Test, HDF5
import DSGE: RootInverseGamma, Exponential
include("../../util.jl")

path = dirname(@__FILE__)

### Model
model = Model990()

### Parameters

# Parameters match parameters, bounds, etc. vectors from reference (ϵ = 1e-4)
para = zeros(82)
bounds = zeros(82, 2)
pshape = zeros(82)
pmean = zeros(82)
pstdd = zeros(82)
trspec = zeros(82, 4)

# keys to skip (used to be fixed_parameters)
fixed_parameters = [:δ, :λ_w, :ϵ_p, :ϵ_w, :g_star]
    
# not all parameters appear in model.parameters
i = 1
for θ in model.parameters
    !isa(θ,AbstractParameter) && continue
    in(θ.key, fixed_parameters) && continue
    
    para[i] = θ.value

    (left, right) = θ.valuebounds
    bounds[i, 1] = left
    bounds[i, 2] = right

    prior = θ.prior.value
    
    if isa(prior, RootInverseGamma)
        pshape[i] = 4
        (ν, τ) = params(prior)
        pmean[i] = τ
        pstdd[i] = ν
    else
        if isa(prior, Distributions.Beta)
            pshape[i] = 1
        elseif isa(prior, Distributions.Gamma)
            pshape[i] = 2
        elseif isa(prior, Distributions.Normal)
            pshape[i] = 3
        end
        pmean[i] = mean(prior)
        pstdd[i] = std(prior)
        
    end

    if θ.transform == Untransformed()
        trspec[i, 1] = 0
    elseif θ.transform == SquareRoot()
        trspec[i, 1] = 1
    elseif  θ.transform == Exponential()
        trspec[i, 1] = 2        
    else
       throw(error("This kind of transform not allowed")) 
    end
        
    (left, right) = θ.transform_parameterization
    trspec[i, 2] = left
    trspec[i, 3] = right
    if θ == model[:modelα_ind]
        trspec[i, 4] = 0
    else
        trspec[i, 4] = 1
    end

    i += 1
end

### Model indices

# Endogenous states
endo = model.endogenous_states
@test length(endo) == 66
@test endo[:E_z] == 60

# Exogenous shocks
exo = model.exogenous_shocks
@test length(exo) == 22
@test exo[:pce_sh] == 16

# Expectation shocks
ex = model.expected_shocks
@test length(ex) == 13
@test ex[:Erk_f_sh] == 13

# Equations
eq = model.equilibrium_conditions
@test length(eq) == 66
@test eq[:eq_Ez] == 60

# Additional states
endo_addl = model.endogenous_states_postgensys
@test length(endo_addl) == 12
@test endo_addl[:y_t1] == 67

# Observables
obs = model.observables
@test length(obs) == 18
@test obs[:tfp] == 12

### Equilibrium conditions
Γ0, Γ1, C, Ψ, Π = eqcond(model)

# Matrices are of expected dimensions
@test size(Γ0) == (66, 66)
@test size(Γ1) == (66, 66)
@test size(C) == (66, 1)
@test size(Ψ) == (66, 22)
@test size(Π) == (66, 13)

# Check output matrices against reference output (ϵ = 1e-4)
h5 = h5open("$path/eqcond.h5")
Γ0_ref = read(h5, "G0")
Γ1_ref = read(h5, "G1")
C_ref  = reshape(read(h5, "C"), 66, 1)
Ψ_ref  = read(h5, "PSI")
Π_ref  = read(h5, "PIE")
close(h5)

@test test_matrix_eq(Γ0_ref, Γ0)
@test test_matrix_eq(Γ1_ref, Γ1)
@test test_matrix_eq(C_ref, C)
@test test_matrix_eq(Ψ_ref, Ψ)
@test test_matrix_eq(Π_ref, Π)


# ### Measurement equation
h5 = h5open("$path/measurement.h5")
ZZ_expected = read(h5, "ZZ")
DD_expected = reshape(read(h5, "DD"), 18, 1)
QQ_expected = read(h5, "QQ")
EE_expected = read(h5, "EE")
MM_expected = read(h5, "MM")
close(h5)

model = Model990()
TTT, RRR, CCC = solve(model)
ZZ, DD, QQ, EE, MM = measurement(model, TTT, RRR, CCC)

@test test_matrix_eq(ZZ_expected, ZZ)
@test test_matrix_eq(DD_expected, DD)
@test test_matrix_eq(QQ_expected, QQ)
@test test_matrix_eq(EE_expected, EE)
@test test_matrix_eq(MM_expected, MM)
