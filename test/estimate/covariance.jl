using Base: Test
using DSGE
using HDF5

include("../util.jl")
path = dirname(@__FILE__)


h5 = h5open("$path/covariance.h5")
hessian = read(h5, "hessian")
U_exp = read(h5, "u")
S_exp = read(h5, "s")
S_inv_exp = read(h5, "s_inv")
cov_exp = read(h5, "invhhm")
close(h5)


S_diag, U = eig(hessian)
S_inv = zeros(82, 82)
for i = 19:82
    S_inv[i, i] = 1/S_diag[i]
end
cov = U*S_inv*U'

@test test_matrix_eq(S_inv_exp, S_inv; verbose=true)
@test test_matrix_eq(cov_exp, cov; verbose=true)
