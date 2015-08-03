using Base: Test
using MATLAB

using DSGE
include("../util.jl")
path = dirname(@__FILE__)



mf = MatFile("$path/svd.mat")
hessian = get_variable(mf, "hessian")
mode = get_variable(mf, "params")
Σ_exp = get_variable(mf, "s")
logdet_exp = get_variable(mf, "sigproplndet")

u_exp = get_variable(mf, "u")
v_exp = get_variable(mf, "v")
σ_exp = get_variable(mf, "sigscale")
close(mf)



# Resulting covariance matrix matches Matlab
d = proposal_distribution(mode, hessian)
@test test_matrix_eq(Σ_exp, d.Σ; noisy=true)
@test_approx_eq logdet_exp d.logdet

# Julia SVD produces a valid decomposition
u, s_diag, v = svd(hessian)
@test test_matrix_eq(hessian, u*diagm(s_diag)*v'; noisy=true)

# Julia SVD with thin=false produces same decomposition because hessian is square
u_thick, s_diag_thick, v_thick = svd(hessian, thin=false)
@test test_matrix_eq(u, u_thick; noisy=true)
@test test_matrix_eq(s_diag, s_diag_thick; noisy=true)
@test test_matrix_eq(v, v_thick; noisy=true)

# u and v matrices not the same as Matlab
@test !test_matrix_eq(u_exp, u; noisy=true)
@test !test_matrix_eq(v_exp, v; noisy=true)

# Also not the same in absolute value
@test !test_matrix_eq(abs(u_exp), abs(u); noisy=true)
@test !test_matrix_eq(abs(v_exp), abs(v); noisy=true)

# Matlab SVD uses LAPACK GESVD
# Julia SVD uses LAPACK GESDD
# What if we call GESVD in Julia?
u_gesvd, s_diag_gesvd, v_gesvd = Base.LinAlg.LAPACK.gesvd!('A', 'A', copy(hessian))
@test !test_matrix_eq(u_exp, u_gesvd; noisy=true)
@test !test_matrix_eq(u, u_gesvd; noisy=true)
@test !test_matrix_eq(abs(u_exp), abs(u_gesvd); noisy=true)
@test !test_matrix_eq(abs(u), abs(u_gesvd); noisy=true)

@test test_matrix_eq([ 0.0 -1.0  0.0],   u_exp[80, 80:82])
@test test_matrix_eq([-1.0  0.0  0.0],       u[80, 80:82])
@test test_matrix_eq([ 0.0  0.0  1.0], u_gesvd[80, 80:82])

# Standard deviation of proposal distribution is cc*σ = cc*u*sqrt(Σ)
jump_exp = σ_exp * ones(82)
jump = d.σ * ones(82)
@test !test_matrix_eq(jump_exp, jump; noisy=true)
@test !test_matrix_eq(abs(jump_exp), abs(jump); noisy=true)
