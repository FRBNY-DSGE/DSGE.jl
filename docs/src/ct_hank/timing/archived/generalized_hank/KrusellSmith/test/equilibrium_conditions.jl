using DSGE, JLD, DataStructures
import Base.Test: @test, @testset

include("../equilibrium_conditions.jl")

# Load required values
params = load("../saved_outputs/test/saved_params.jld", "params")
grids = load("../saved_outputs/test/saved_params.jld", "grids")
init_params = load("../saved_outputs/test/saved_params.jld", "init_params")
approx_params = load("../saved_outputs/test/saved_params.jld", "approx_params")
varsSS = load("../saved_outputs/test/varsSS.jld", "varsSS")

# Load true canonical matrices
g0 = load("../saved_outputs/test/canonical.jld", "g0")
g1 = load("../saved_outputs/test/canonical.jld", "g1")
pi = load("../saved_outputs/test/canonical.jld", "pi")
psi = load("../saved_outputs/test/canonical.jld", "psi")
c = load("../saved_outputs/test/canonical.jld", "c")

# Test
canonical_form = equilibrium_conditions(varsSS, grids, params, approx_params)
@testset "Derivatives of equilibrium conditions" begin
    @test @test_matrix_approx_eq full(g0) canonical_form[:Γ0]
    @test @test_matrix_approx_eq full(g1) canonical_form[:Γ1]
    @test vec(full(pi)) ≈ vec(full(canonical_form[:Π]))
    @test vec(full(psi)) ≈ vec(full(canonical_form[:Ψ]))
    @test vec(full(c)) ≈ vec(full(canonical_form[:C]))
end