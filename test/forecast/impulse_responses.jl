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

# Add parallel workers
ndraws = 2
my_procs = addprocs(ndraws)
@everywhere using DSGE

systems = distribute(systems; procs = my_procs, dist = [ndraws])

# Run to compile before timing
_, _, _ = DSGE.impulse_responses(m, systems; procs = my_procs)

# Run and time
@time states, obs, pseudo = DSGE.impulse_responses(m, systems; procs = my_procs)

# Compare to expected output
exp_states, exp_obs, exp_pseudo = jldopen("$path/../reference/impulse_responses_out.jld", "r") do file
    read(file, "exp_states"), read(file, "exp_obs"), read(file, "exp_pseudo")
end

@test_matrix_approx_eq exp_states convert(Array, states)
@test_matrix_approx_eq exp_obs    convert(Array, obs)
@test_matrix_approx_eq exp_pseudo convert(Array, pseudo)

# Remove parallel workers
rmprocs(my_procs)

nothing
