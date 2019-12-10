using DSGE, JLD, DataStructures
import Base.Test: @test, @testset

include("../../ReductionData_struct.jl")
include("../../reduction.jl")

# load parameters
red_params = load("../saved_outputs/test/saved_params.jld", "red_params")
grids = load("../saved_outputs/test/saved_params.jld", "grids")

# Load pre-reduced canonical form matrices
g1 = load("../saved_outputs/test/canonical.jld", "g1")
g0 = load("../saved_outputs/test/canonical.jld", "g0")
psi = load("../saved_outputs/test/canonical.jld", "psi")
pi = load("../saved_outputs/test/canonical.jld", "pi")
c = load("../saved_outputs/test/canonical.jld", "c")

# Load reduced canonical form matrices
g1_kry = load("../saved_outputs/test/canonical_kry.jld", "g1")
g0_kry = load("../saved_outputs/test/canonical_kry.jld", "g0")
psi_kry = load("../saved_outputs/test/canonical_kry.jld", "psi")
pi_kry = load("../saved_outputs/test/canonical_kry.jld", "pi")
c_kry = load("../saved_outputs/test/canonical_kry.jld", "c")

# Test
Γ0, Γ1, Ψ, Π, C, ~, ~ = krylov_reduction(full(g0), full(g1), full(psi), full(pi), full(c), grids, red_params)
@testset "Krylov Reduction" begin
    @test @test_matrix_approx_eq full(g1_kry) Γ1
    @test @test_matrix_approx_eq full(g0_kry) Γ0
    @test vec(full(psi_kry)) ≈ vec(full(Ψ)) # allow machine error
    @test vec(full(pi_kry)) ≈ vec(full(Π))
    @test vec(full(c_kry)) ≈ vec(full(C))
end
