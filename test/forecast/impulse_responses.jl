using DSGE, JLD

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

@test_matrix_approx_eq exp_states states
@test_matrix_approx_eq exp_obs    obs
@test_matrix_approx_eq exp_pseudo pseudo

nothing
