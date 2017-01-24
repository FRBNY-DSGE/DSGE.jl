using DSGE, JLD, DistributedArrays
include("../util.jl")

path = dirname(@__FILE__)

# Set up arguments
custom_settings = Dict{Symbol, Setting}(
    :n_anticipated_shocks => Setting(:n_anticipated_shocks, 6),
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2016-Q1")),
    :forecast_horizons    => Setting(:forecast_horizons, 1),
    :use_parallel_workers => Setting(:use_parallel_workers, true),
    :forecast_kill_shocks => Setting(:forecast_kill_shocks, true))
m = Model990(custom_settings = custom_settings, testing = true)

systems, z0s, kals = jldopen("$path/../reference/forecast_args.jld","r") do file
    read(file, "systems"), read(file, "states"), read(file, "kals")
end

# Add parallel workers
ndraws = length(systems) # 2
my_procs = addprocs(ndraws)
@everywhere using DSGE

# Distribute systems and states
systems = distribute(systems; procs = my_procs, dist = [ndraws])
z0s     = distribute(z0s;     procs = my_procs, dist = [ndraws])
kals    = distribute(kals;    procs = my_procs, dist = [ndraws])

# Run to compile before timing
states, obs, pseudo, shocks = DSGE.forecast(m, systems, z0s; procs = my_procs)
states, obs, pseudo, shocks = DSGE.forecast(m, systems, kals; procs = my_procs)

# Read expected output
exp_states, exp_obs, exp_pseudo, exp_shocks = jldopen("$path/../reference/forecast_out.jld", "r") do file
    read(file, "exp_states"), read(file, "exp_obs"), read(file, "exp_pseudo"), read(file, "exp_shocks")
end

# Run forecast without supplying shocks
@time states, obs, pseudo, shocks = DSGE.forecast(m, systems, z0s; procs = my_procs)

@test_matrix_approx_eq exp_states convert(Array, states)
@test_matrix_approx_eq exp_obs    convert(Array, obs)
@test_matrix_approx_eq exp_pseudo convert(Array, pseudo)
@test_matrix_approx_eq exp_shocks convert(Array, shocks)

# Run forecast, supplying shocks
@time states, obs, pseudo, shocks = DSGE.forecast(m, systems, z0s; procs = my_procs, shocks = shocks)

@test_matrix_approx_eq exp_states convert(Array, states)
@test_matrix_approx_eq exp_obs    convert(Array, obs)
@test_matrix_approx_eq exp_pseudo convert(Array, pseudo)
@test_matrix_approx_eq exp_shocks convert(Array, shocks)

# Run forecast, drawing z0s
m <= Setting(:forecast_draw_z0, true)
@time states, obs, pseudo, shocks = DSGE.forecast(m, systems, kals; procs = my_procs)

# Remove parallel workers
rmprocs(my_procs)

nothing
