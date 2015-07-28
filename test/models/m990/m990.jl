using Distributions
using DSGE: DistributionsExt
path = dirname(@__FILE__)

model = Model990()
@test isa(model, Model990)

# Parameters990 object creation
Θ = Parameters990(model_specifications(Model990))
@test isa(Θ, Parameters990)

# Parameters match para, bounds, etc. vectors from Matlab (ε = 1e-4)
para = zeros(82, 1)
bounds = zeros(82, 2)
pshape = zeros(82, 1)
pmean = zeros(82, 1)
pstdd = zeros(82, 1)
trspec = zeros(82, 4)

# not all Params appear in para vector
i = 1
for φ = Θ
    para[i] = φ.value

    (left, right) = φ.bounds
    bounds[i, 1] = left
    bounds[i, 2] = right

    if isa(φ.priordist, DistributionsExt.RootInverseGamma)
        pshape[i] = 4
        (ν, τ) = params(φ.priordist)
        pmean[i] = τ
        pstdd[i] = ν
    else
        if isa(φ.priordist, Distributions.Beta)
            pshape[i] = 1
        elseif isa(φ.priordist, Distributions.Gamma)
            pshape[i] = 2
        elseif isa(φ.priordist, Distributions.Normal)
            pshape[i] = 3
        end
        pmean[i] = mean(φ.priordist)
        pstdd[i] = std(φ.priordist)
    end

    trspec[i, 1] = φ.transformtype
    (left, right) = φ.transformbounds
    trspec[i, 2] = left
    trspec[i, 3] = right
    if φ == Θ.modelalp_ind
        trspec[i, 4] = 0
    else
        trspec[i, 4] = 1
    end

    i += 1
end


para_matlab   = readcsv(joinpath(path,"parameters/para.csv"))
bounds_matlab = readcsv(joinpath(path,"parameters/bounds.csv"))
pshape_matlab = readcsv(joinpath(path,"parameters/pshape.csv"))
pmean_matlab  = readcsv(joinpath(path,"parameters/pmean.csv"))
pstdd_matlab  = readcsv(joinpath(path,"parameters/pstdd.csv"))
trspec_matlab = readcsv(joinpath(path,"parameters/trspec.csv"))

@test test_matrix_eq(para_matlab, para)
@test test_matrix_eq(bounds_matlab, bounds)
@test test_matrix_eq(pshape_matlab, pshape)
@test test_matrix_eq(pmean_matlab, pmean)
@test test_matrix_eq(pstdd_matlab, pstdd)
@test test_matrix_eq(trspec_matlab, trspec)

# ModelInds object creation
I = ModelInds(model_specifications(Model990))
@test isa(I, ModelInds)

# Endogenous states
endo = I.endostates
@test length(endo) == 66
@test endo["E_z"] == 60

# Exogenous shocks
exo = I.exoshocks
@test length(exo) == 22
@test exo["pce_sh"] == 16

# Expectation shocks
ex = I.expshocks
@test length(ex) == 13
@test ex["Erk_f_sh"] == 13

# Equations
eq = I.eqconds
@test length(eq) == 66
@test eq["eq_Ez"] == 60

# Additional states
endo_addl = I.endostates_postgensys
@test length(endo_addl) == 12
@test endo_addl["y_t1"] == 67

# Observables
obs = I.observables
@test length(obs) == 18
@test obs["tfp"] == 12

model = Model990()
Γ0, Γ1, C, Ψ, Π = eqcond(model)

# Matrices are of expected dimensions
@test size(Γ0) == (66, 66)
@test size(Γ1) == (66, 66)
@test size(C) == (66, 1)
@test size(Ψ) == (66, 22)
@test size(Π) == (66, 13)

# Check output matrices against Matlab output (ε = 1e-4)
Γ0_matlab = readcsv(joinpath(path,"eqcond/Γ0.csv"))
Γ1_matlab = readcsv(joinpath(path,"eqcond/Γ1.csv"))
C_matlab  = readcsv(joinpath(path,"eqcond/C.csv"))
Ψ_matlab  = readcsv(joinpath(path,"eqcond/PSI.csv"))
Π_matlab  = readcsv(joinpath(path,"eqcond/PIE.csv"))

@test test_matrix_eq(Γ0_matlab, Γ0)
@test test_matrix_eq(Γ1_matlab, Γ1)
@test test_matrix_eq(C_matlab, C)
@test test_matrix_eq(Ψ_matlab, Ψ)
@test test_matrix_eq(Π_matlab, Π)

# MATLAB Code

using MATLAB

mf = MatFile("$path/test_measurement.mat")
ZZ_expected = get_variable(mf, "ZZ")
DD_expected = reshape(get_variable(mf, "DD"), 18, 1)
QQ_expected = get_variable(mf, "QQ")
EE_expected = get_variable(mf, "EE")
MM_expected = get_variable(mf, "MM")
close(mf)

model = Model990()
TTT, RRR, CCC = solve(model)
ZZ, DD, QQ, EE, MM = measurement(model, TTT, RRR, CCC)

@test test_matrix_eq(ZZ_expected, ZZ)
@test test_matrix_eq(DD_expected, DD)
@test test_matrix_eq(QQ_expected, QQ)
@test test_matrix_eq(EE_expected, EE)
@test test_matrix_eq(MM_expected, MM)
