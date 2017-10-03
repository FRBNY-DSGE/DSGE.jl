using HDF5, Base.Test
using ClusterManagers

# Read in an input matrix and reference eigenvalues/eigenvectors
h5 = h5open(joinpath(dirname(@__FILE__), "..", "reference", "eig.h5"),"r")
matrix = read(h5,"matrix_in")
ref_U = read(h5,"U")
ref_S_diag = read(h5, "S_diag")
close(h5)

# Eigendecomposition with only 1 worker
my_proc = addprocs_sge(1)
S_diag, U = eig(matrix)

# Test for equality
@test_matrix_approx_eq_eps ref_S_diag S_diag 1e-6 1e-2
@test_matrix_approx_eq_eps ref_U U 1e-6 1e-2
# TODO calibrate above for when use_parallel_workers ≡ true

rmprocs(my_proc)

nothing
