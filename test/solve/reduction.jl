using DSGE, JLD, DataStructures
import Base.Test: @test, @testset
import DSGE: @test_matrix_approx_eq

m = KrusellSmith()

# Load pre-reduced canonical form matrices
g1 = load("../../test_outputs/solve/canonical.jld", "g1")
g0 = load("../../test_outputs/solve/canonical.jld", "g0")
psi = load("../../test_outputs/solve/canonical.jld", "psi")
pi = load("../../test_outputs/solve/canonical.jld", "pi")
c = load("../../test_outputs/solve/canonical.jld", "c")

# Load reduced canonical form matrices
g1_kry = load("../../test_outputs/solve/canonical_kry.jld", "g1")
g0_kry = load("../../test_outputs/solve/canonical_kry.jld", "g0")
psi_kry = load("../../test_outputs/solve/canonical_kry.jld", "psi")
pi_kry = load("../../test_outputs/solve/canonical_kry.jld", "pi")
c_kry = load("../../test_outputs/solve/canonical_kry.jld", "c")

# Load Krylov + spline reduced canonical form matrices
g1_spl = load("../../test_outputs/solve/canonical_spl.jld", "g1")
g0_spl = load("../../test_outputs/solve/canonical_spl.jld", "g0")
psi_spl = load("../../test_outputs/solve/canonical_spl.jld", "psi")
pi_spl = load("../../test_outputs/solve/canonical_spl.jld", "pi")
c_spl = load("../../test_outputs/solve/canonical_spl.jld", "c")

# Update model
EHANK.update!(m.settings[:I], Setting(:I, 100))
EHANK.update!(m.settings[:amin], Setting(:amin, 0.))
EHANK.update!(m.settings[:amax], Setting(:amax, 100.))
EHANK.update!(m.settings[:a], Setting(:a, collect(linspace(0., 100., 100))))
EHANK.update!(m.settings[:da], Setting(:da, 100./(100 - 1)))
EHANK.update!(m.settings[:aaa], Setting(:aaa, reshape(get_setting(m, :aa), 2 * 100, 1)))
EHANK.update!(m.settings[:J], Setting(:J, 2))
EHANK.update!(m.settings[:z], Setting(:z, [0., 1.]))
EHANK.update!(m.settings[:dz], Setting(:dz, 1.))
EHANK.update!(m.settings[:zz], Setting(:zz, ones(100, 1) * get_setting(m, :z)'))
EHANK.update!(m.settings[:zzz], Setting(:zzz, reshape(get_setting(m, :zz), 2 * 100, 1)))
EHANK.update!(m.settings[:n_jump_vars], Setting(:n_jump_vars, 200))
EHANK.update!(m.settings[:n_state_vars], Setting(:n_state_vars, 199 + 1))
EHANK.update!(m.settings[:n_state_vars_unreduce], Setting(:n_state_vars_unreduce, 0))
EHANK.update!(m.settings[:reduce_state_vars], Setting(:reduce_state_vars, true))
EHANK.update!(m.settings[:reduce_v], Setting(:reduce_v, true))
EHANK.update!(m.settings[:krylov_dim], Setting(:krylov_dim, 5))
EHANK.update!(m.settings[:n_knots], Setting(:n_knots, 12))
EHANK.update!(m.settings[:c_power], Setting(:c_power, 7))
EHANK.update!(m.settings[:knots_dict], Setting(:knots_dict, Dict(1 => collect(linspace(get_setting(m, :amin), get_setting(m, :amax), 12 - 1)))))
EHANK.update!(m.settings[:spline_grid], Setting(:spline_grid, get_setting(m, :a)))
EHANK.update!(m.settings[:n_prior], Setting(:n_prior, 1))
EHANK.update!(m.settings[:n_post], Setting(:n_post, 2))
EHANK.update!(m.settings[:F], Setting(:F, identity))

# Test non-sparse methods
Γ0, Γ1, Ψ, Π, C, ~, ~ = krylov_reduction(m, full(g0), full(g1), full(psi), full(pi), full(c))
@testset "Krylov Reduction (non-sparse)" begin
    @test @test_matrix_approx_eq full(g1_kry) Γ1
    @test @test_matrix_approx_eq full(g0_kry) Γ0
    @test vec(full(psi_kry)) ≈ vec(full(Ψ)) # allow machine error
    @test vec(full(pi_kry)) ≈ vec(full(Π))
    @test vec(full(c_kry)) ≈ vec(full(C))
end

Γ0, Γ1, Ψ, Π, C, ~, ~ = valuef_reduction(m, full(g0_kry), full(g1_kry), full(psi_kry), full(pi_kry), full(c_kry))
@testset "Value Function Reduction via Spline Basis (non-sparse)" begin
    @test @test_matrix_approx_eq full(g1_spl) Γ1
    @test @test_matrix_approx_eq full(g0_spl) Γ0
    @test vec(full(psi_spl)) ≈ vec(full(Ψ))
    @test vec(full(pi_spl)) ≈ vec(full(Π))
    @test vec(full(c_spl)) ≈ vec(full(C))
end

# Test sparse methods
Γ0, Γ1, Ψ, Π, C, ~, ~ = krylov_reduction(m, sparse(g0), sparse(g1), sparse(psi), sparse(pi), sparse(c))
@testset "Krylov Reduction (sparse)" begin
    @test typeof(Γ1) == SparseMatrixCSC{Float64,Int64}
    @test typeof(Γ0) == SparseMatrixCSC{Float64,Int64}
    @test typeof(Ψ) == SparseMatrixCSC{Float64,Int64}
    @test typeof(Π) == SparseMatrixCSC{Float64,Int64}
    @test typeof(C) == SparseMatrixCSC{Float64,Int64}
    @test @test_matrix_approx_eq full(g1_kry) full(Γ1)
    @test @test_matrix_approx_eq full(g0_kry) full(Γ0)
    @test vec(full(psi_kry)) ≈ vec(full(Ψ)) # allow machine error
    @test vec(full(pi_kry)) ≈ vec(full(Π))
    @test vec(full(c_kry)) ≈ vec(full(C))
end

Γ0, Γ1, Ψ, Π, C, ~, ~ = valuef_reduction(m, sparse(g0_kry), sparse(g1_kry), sparse(psi_kry),
                                         sparse(pi_kry), sparse(c_kry))
@testset "Value Function Reduction via Spline Basis (sparse)" begin
    @test typeof(Γ1) == SparseMatrixCSC{Float64,Int64}
    @test typeof(Γ0) == SparseMatrixCSC{Float64,Int64}
    @test typeof(Ψ) == SparseMatrixCSC{Float64,Int64}
    @test typeof(Π) == SparseMatrixCSC{Float64,Int64}
    @test typeof(C) == SparseMatrixCSC{Float64,Int64}
    @test @test_matrix_approx_eq full(g1_spl) full(Γ1)
    @test @test_matrix_approx_eq full(g0_spl) full(Γ0)
    @test vec(full(psi_spl)) ≈ vec(full(Ψ))
    @test vec(full(pi_spl)) ≈ vec(full(Π))
    @test vec(full(c_spl)) ≈ vec(full(C))
end
