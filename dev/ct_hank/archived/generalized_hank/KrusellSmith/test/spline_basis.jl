using DSGE, JLD, BasisMatrices
import Base.Test: @test, @testset

include("../../ReductionData_struct.jl")
include("../../reduction.jl")

# Load parameters
grids = load("../saved_outputs/test/saved_params.jld", "grids")
red_params = load("../saved_outputs/test/saved_params.jld", "red_params")

# Load Krylov reduced canonical form matrices
g1_kry = load("../saved_outputs/test/canonical_kry.jld", "g1")
g0_kry = load("../saved_outputs/test/canonical_kry.jld", "g0")
psi_kry = load("../saved_outputs/test/canonical_kry.jld", "psi")
pi_kry = load("../saved_outputs/test/canonical_kry.jld", "pi")
c_kry = load("../saved_outputs/test/canonical_kry.jld", "c")

# Load Krylov + spline reduced canonical form matrices
g1_spl = load("../saved_outputs/test/canonical_spl.jld", "g1")
g0_spl = load("../saved_outputs/test/canonical_spl.jld", "g0")
psi_spl = load("../saved_outputs/test/canonical_spl.jld", "psi")
pi_spl = load("../saved_outputs/test/canonical_spl.jld", "pi")
c_spl = load("../saved_outputs/test/canonical_spl.jld", "c")

# Test
Γ0, Γ1, Ψ, Π, C, ~, ~ = valuef_reduction(grids, red_params, full(g0_kry), full(g1_kry), full(psi_kry), full(pi_kry), full(c_kry))
@testset "Value Function Reduction via Spline Basis" begin
    @test @test_matrix_approx_eq full(g1_spl) Γ1
    @test @test_matrix_approx_eq full(g0_spl) Γ0
    @test vec(full(psi_spl)) ≈ vec(Ψ)
    @test vec(full(pi_spl)) ≈ vec(Π)
    @test vec(full(c_spl)) ≈ vec(C)
end
