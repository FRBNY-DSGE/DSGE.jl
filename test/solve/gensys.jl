using Base: Test
using HDF5
using DSGE

include("../util.jl")
path = dirname(@__FILE__)

m = AnSchorfheide()
Γ0, Γ1, C, Ψ, Π = eqcond(m)
stake = 1 + 1e-6
G1, C, impact, fmat, fwt, ywt, gev, eu, loose = gensys(Γ0, Γ1, C, Ψ, Π, stake)

file = "$path/../reference/gensys.h5"
G1_exp = h5read(file, "G1_gensys")
C_exp = h5read(file,"C_gensys")
impact_exp = h5read(file, "impact")
eu_exp = h5read(file, "eu")

@test_matrix_approx_eq G1_exp G1

@test_matrix_approx_eq C_exp C

@test_matrix_approx_eq impact_exp impact

@test isequal(eu_exp, eu)
