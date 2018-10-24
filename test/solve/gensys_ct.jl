using DSGE, JLD, DataStructures
import Test: @test, @testset
import DSGE: @test_matrix_approx_eq

# Load pre-input files
out = load("../../test_outputs/solve/splines_outputs_mat.jld")
g0 = full(out["g0"]); g1 = full(out["g1"]); psi = full(out["psi"])
Pi = full(out["pi"]); c = full(out["constant"])
out = load("../../test_outputs/solve/splines_outputs.jld")
Γ0 = full(out["Γ0"]); Γ1 = full(out["Γ1"]); Ψ = full(out["Ψ"])
Π = full(out["Π"]); C_val = full(out["C"])

# Load checking files
out = load("../../test_outputs/solve/post_gensys_mat.jld")
G1_mat = full(out["G1"]); impact_mat = full(out["impact"])
out = load("../../test_outputs/solve/post_gensys.jld")
G1 = full(out["G1"]); impact = full(out["impact"]);

# Check non-in-place first
T_mat, ~, R_mat, ~, ~, eu_mat = gensysct(g1, c, psi, Pi, complex_decomposition = false)
T, ~, R, ~, ~, ~, ~, eu = gensysct(Γ0, Γ1, C_val, Ψ, Π, complex_decomposition = false)
@testset "Non-in-place gensysct" begin
    @test @test_matrix_approx_eq T_mat G1_mat
    @test @test_matrix_approx_eq R_mat impact_mat
    @test eu_mat == [1; 1]
    @test @test_matrix_approx_eq T G1
    @test @test_matrix_approx_eq R impact
    @test eu == [-1; 1]
end

# Check in place gensysct
T_mat, ~, R_mat, ~, ~, eu_mat = gensysct!(g1, c, psi, Pi, complex_decomposition = false)
T, ~, R, ~, ~, ~, ~, eu = gensysct!(Γ0, Γ1, C_val, Ψ, Π, complex_decomposition = false)
@testset "In-place gensysct" begin
    @test @test_matrix_approx_eq T_mat G1_mat
    @test @test_matrix_approx_eq R_mat impact_mat
    @test eu_mat == [1; 1]
    @test @test_matrix_approx_eq T G1
    @test @test_matrix_approx_eq R impact
    @test eu == [-1; 1]
end