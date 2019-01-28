using DSGE, JLD
using Base.Test

path = dirname(@__FILE__)

# Set up arguments
m = AnSchorfheide(testing = true)

system = jldopen("$path/../reference/forecast_args.jld","r") do file
    read(file, "system")
end

# Run impulse responses
states, obs, pseudo = impulse_responses(m, system)

# Compare to expected output
exp_states, exp_obs, exp_pseudo =
    jldopen("$path/../reference/impulse_responses_out.jld", "r") do file
        read(file, "exp_states"), read(file, "exp_obs"), read(file, "exp_pseudo")
    end

@testset "Compare irfs to expected output" begin
    @test @test_matrix_approx_eq exp_states states
    @test @test_matrix_approx_eq exp_obs    obs
    @test @test_matrix_approx_eq exp_pseudo pseudo
end

# Test impulse response methods for:
# (1) Choosing a set of shocks and their shock values
# (2) Choosing a shock, a state, and the value the state should hit in period 1 based on
# an impulse from that shock at time 1.
# (3) Same as (2) but with an observable.

# Setup
endo   = m.endogenous_states
obsvs  = m.observables
exo    = m.exogenous_shocks
horizon = 10

# Testing both for exact output match and for certain values to match to be
# explicit about the expected outcome of the impulse response

# Shock
shock_names  = [:z_sh, :g_sh]
shock_values = [1., 1.]
states, obs, pseudo = impulse_responses(m, system, horizon, shock_names, shock_values)

# Explicit checks
@testset "Compare shock-set irfs to the states they affect 1-for-1" begin
    @test states[endo[:z_t], 1, exo[:z_sh]] ≈ 1.
    @test states[endo[:y_t], 1, exo[:g_sh]] ≈ 1.
    @test states[endo[:g_t], 1, exo[:g_sh]] ≈ 1.
end

exp_states_shockset, exp_obs_shockset, exp_pseudo_shockset =
    jldopen("$path/../reference/impulse_responses_out.jld", "r") do file
        read(file, "exp_states_shockset"), read(file, "exp_obs_shockset"), read(file, "exp_pseudo_shockset")
    end

@testset "Checking shockset irfs" begin
    @test states ≈ exp_states_shockset
    @test obs    ≈ exp_obs_shockset
    @test pseudo ≈ exp_pseudo_shockset
end

# State
shock_name = :z_sh
state_name   = :y_t
state_value  = 1.
states, obs, pseudo = impulse_responses(m, system, horizon, shock_name, state_name, state_value)

# Explicit checks
@testset "Check state irf initial period is exactly the specified state value" begin
    @test states[endo[:y_t], 1, exo[:z_sh]] ≈ 1.
end

exp_states_shockstates, exp_obs_shockstates, exp_pseudo_shockstates =
    jldopen("$path/../reference/impulse_responses_out.jld", "r") do file
        read(file, "exp_states_shockstate"), read(file, "exp_obs_shockstate"), read(file, "exp_pseudo_shockstate")
    end

@testset "Checking shockstate irfs" begin
    @test states ≈ exp_states_shockstates
    @test obs    ≈ exp_obs_shockstates
    @test pseudo ≈ exp_pseudo_shockstates
end

# Observable
shock_name = :z_sh
obs_name   = :obs_gdp
obs_value  = 1.
states, obs, pseudo = impulse_responses(m, system, horizon, shock_name, obs_name, obs_value)

# Explicit checks
@testset "Check obs irf initial period is exactly the specified obs value" begin
    @test obs[obsvs[:obs_gdp], 1, exo[:z_sh]] ≈ 1.
end

exp_states_shockobs, exp_obs_shockobs, exp_pseudo_shockobs =
    jldopen("$path/../reference/impulse_responses_out.jld", "r") do file
        read(file, "exp_states_shockobs"), read(file, "exp_obs_shockobs"), read(file, "exp_pseudo_shockobs")
    end

@testset "Checking shockobs irfs" begin
    @test states ≈ exp_states_shockobs
    @test obs    ≈ exp_obs_shockobs
    @test pseudo ≈ exp_pseudo_shockobs
end

nothing
