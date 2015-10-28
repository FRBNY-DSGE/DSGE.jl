using Base.Test, HDF5
include("../util.jl")

# Read in an input matrix and reference eigenvalues/eigenvectors
h5 = h5open(joinpath(dirname(@__FILE__), "..", "reference", "eig.h5"),"r")
matrix = read(h5,"matrix_in")
ref_U = read(h5,"U")
ref_S_diag = read(h5, "S_diag")
close(h5)

# Eigendecomposition with only 1 worker
addprocs(1) 
S_diag, U = eig(matrix) 

# Test for equality
@test test_matrix_eq(ref_S_diag, S_diag, ϵ=1e-9, ϵ_pct=1e-9)
@test test_matrix_eq(ref_U, U, ϵ=1e-9, ϵ_pct=1e-9)

