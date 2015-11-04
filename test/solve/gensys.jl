using Base: Test, LinAlg
using HDF5

using DSGE
include("../util.jl")
path = dirname(@__FILE__)



model = Model990()
Γ0, Γ1, C, Ψ, Π = eqcond(model)
stake = 1 + 1e-6
G1, C, impact, fmat, fwt, ywt, gev, eu, loose = gensys(Γ0, Γ1, C, Ψ, Π, stake)

h5 = h5open("$path/../reference/gensys.h5")
G1_exp = read(h5, "G1_gensys")
C_exp = reshape(read(h5, "C_gensys"), 66, 1)
impact_exp = read(h5, "impact")
eu_exp = read(h5, "eu")
close(h5)

@test_matrix_approx_eq G1_exp G1

@test_matrix_approx_eq C_exp C

@test_matrix_approx_eq impact_exp impact

@test isequal(eu_exp, eu)
