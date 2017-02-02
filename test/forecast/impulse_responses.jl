using DSGE, JLD
include("../util.jl")

path = dirname(@__FILE__)

# Set up arguments
custom_settings = Dict{Symbol, Setting}(
    :n_anticipated_shocks => Setting(:n_anticipated_shocks, 6),
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")),
    :use_parallel_workers => Setting(:use_parallel_workers, true))
m = Model990(custom_settings = custom_settings, testing = true)

system = jldopen("$path/../reference/forecast_args.jld","r") do file
    read(file, "system")
end

# Run impulse responses
states, obs, pseudo = impulse_responses(m, system)

# Compare to expected output
exp_states, exp_obs, exp_pseudo = jldopen("$path/../reference/impulse_responses_out.jld", "r") do file
    read(file, "exp_states"), read(file, "exp_obs"), read(file, "exp_pseudo")
end

@test_matrix_approx_eq squeeze(exp_states[1, :, :, :], 1) states
@test_matrix_approx_eq squeeze(exp_obs[1, :, :, :],    1) obs
@test_matrix_approx_eq squeeze(exp_pseudo[1, :, :, :], 1) pseudo


nothing
