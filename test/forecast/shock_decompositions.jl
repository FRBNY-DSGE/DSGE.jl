using DSGE, JLD, DistributedArrays
include("../util.jl")

path = dirname(@__FILE__)

# Set up arguments
custom_settings = Dict{Symbol, Setting}(
    :n_anticipated_shocks => Setting(:n_anticipated_shocks, 6),
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")),
    :use_parallel_workers => Setting(:use_parallel_workers, true))
m = Model990(custom_settings = custom_settings, testing = true)

systems, histshocks = jldopen("$path/../reference/forecast_args.jld","r") do file
    read(file, "systems"), read(file, "histshocks")
end
ndraws     = length(systems) # 2
systems    = distribute(systems;    procs = [myid()])
histshocks = distribute(histshocks; procs = [myid()])

# Run to compile before timing
shock_decompositions(m, systems, histshocks)

# Read expected output
exp_states, exp_obs, exp_pseudo = jldopen("$path/../reference/shock_decompositions_out.jld", "r") do file
    read(file, "exp_states"), read(file, "exp_obs"), read(file, "exp_pseudo")
end

# With shockdec_startdate not null
@time states, obs, pseudo = shock_decompositions(m, systems, histshocks)

@test_matrix_approx_eq exp_states[:startdate] convert(Array, states)
@test_matrix_approx_eq exp_obs[:startdate]    convert(Array, obs)
@test_matrix_approx_eq exp_pseudo[:startdate] convert(Array, pseudo)

# With shockdec_startdate null
m <= Setting(:shockdec_startdate, Nullable{Date}())
@time states, obs, pseudo = shock_decompositions(m, systems, histshocks)

@test_matrix_approx_eq exp_states[:no_startdate] convert(Array, states)
@test_matrix_approx_eq exp_obs[:no_startdate]    convert(Array, obs)
@test_matrix_approx_eq exp_pseudo[:no_startdate] convert(Array, pseudo)

nothing
