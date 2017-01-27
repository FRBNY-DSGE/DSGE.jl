using DSGE, JLD, DistributedArrays
include("../util.jl")

path = dirname(@__FILE__)

# Set up arguments
custom_settings = Dict{Symbol, Setting}(
    :n_anticipated_shocks => Setting(:n_anticipated_shocks, 6),
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")),
    :use_parallel_workers => Setting(:use_parallel_workers, true))
m = Model990(custom_settings = custom_settings, testing = true)

systems = jldopen("$path/../reference/forecast_args.jld","r") do file
    read(file, "systems")
end
ndraws  = length(systems) # 2
systems = distribute(systems; procs = [myid()])

# Run to compile before timing
impulse_responses(m, systems)

# Run and time
@time states, obs, pseudo = impulse_responses(m, systems)

# Compare to expected output
exp_states, exp_obs, exp_pseudo = jldopen("$path/../reference/impulse_responses_out.jld", "r") do file
    read(file, "exp_states"), read(file, "exp_obs"), read(file, "exp_pseudo")
end

@test_matrix_approx_eq exp_states convert(Array, states)
@test_matrix_approx_eq exp_obs    convert(Array, obs)
@test_matrix_approx_eq exp_pseudo convert(Array, pseudo)


nothing
