using DSGE, JLD, DataStructures
import Base.Test: @test, @testset

include("../compute_steady_state.jl")
include("../../helpers.jl")

# Load parameters
params = load("../saved_outputs/test/saved_params.jld", "params")
grids = load("../saved_outputs/test/saved_params.jld", "grids")
init_params = load("../saved_outputs/test/saved_params.jld", "init_params")
approx_params = load("../saved_outputs/test/saved_params.jld", "approx_params")

# Load true varsSS
vars_SS_mat = load("../saved_outputs/test/varsSS.jld")
vars_SS_mat = vars_SS_mat["varsSS"]

# Test
varsSS = compute_steady_state(grids, params, init_params, approx_params)
@testset "Steady State" begin
    @test @test_matrix_approx_eq vec(vars_SS_mat[:VSS]) vec(varsSS[:VSS])
    @test @test_matrix_approx_eq vec(vars_SS_mat[:ggSS]) vec(varsSS[:ggSS])
    @test vars_SS_mat[:KSS] ≈ varsSS[:KSS] # allow for machine error
    @test vars_SS_mat[:rSS] ≈ varsSS[:rSS]
    @test vars_SS_mat[:wSS] ≈ varsSS[:wSS]
    @test vars_SS_mat[:YSS] ≈ varsSS[:YSS]
    @test vars_SS_mat[:CSS] ≈ varsSS[:CSS]
    @test vars_SS_mat[:ISS] ≈ varsSS[:ISS]
end
