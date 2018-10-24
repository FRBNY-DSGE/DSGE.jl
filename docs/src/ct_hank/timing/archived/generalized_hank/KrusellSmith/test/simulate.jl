using DSGE, JLD
import Base.Test: @test, @testset

include("../../simulate.jl")

# Load values for computing sim_states
periods = load("../saved_outputs/test/sim_states.jld", "periods")
steps = load("../saved_outputs/test/sim_states.jld", "steps")
agg_shock = load("../saved_outputs/test/sim_states.jld", "agg_shock")
inverse_basis = load("../saved_outputs/test/sim_states.jld", "inverse_basis")
method = load("../saved_outputs/test/sim_states.jld", "method")
to_extract = load("../saved_outputs/test/sim_states.jld", "to_extract")
G1 = load("../saved_outputs/test/sim_states.jld", "G1")
impact = load("../saved_outputs/test/sim_states.jld", "impact")

# Load true sim_states
saved_sim_states = load("../saved_outputs/test/sim_states.jld", "sim_states")

# Test. If sim_states is correct, then IRFs are, assuming compute_steady_state went correctly
sim_states = simulate(full(G1), full(impact), periods, steps, agg_shock;
                      method = method, transformation = full(inverse_basis), to_extract = to_extract)
@testset "Simulating Impulse Responses" begin
    @test @test_matrix_approx_eq sim_states saved_sim_states
end
