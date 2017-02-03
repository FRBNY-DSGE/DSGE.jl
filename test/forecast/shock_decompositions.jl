using DSGE, JLD
include("../util.jl")

path = dirname(@__FILE__)

# Set up arguments
custom_settings = Dict{Symbol, Setting}(
    :n_anticipated_shocks => Setting(:n_anticipated_shocks, 6),
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")),
    :use_parallel_workers => Setting(:use_parallel_workers, true))
m = Model990(custom_settings = custom_settings, testing = true)

system, histshocks = jldopen("$path/../reference/forecast_args.jld","r") do file
    read(file, "system"), read(file, "histshocks")
end

# Read expected output
exp_states, exp_obs, exp_pseudo = jldopen("$path/../reference/shock_decompositions_out.jld", "r") do file
    read(file, "exp_states"), read(file, "exp_obs"), read(file, "exp_pseudo")
end

# With shockdec_startdate not null
states, obs, pseudo = shock_decompositions(m, system, histshocks)

# @test_matrix_approx_eq squeeze(exp_states[:startdate][1, :, :, :], 1) states
# @test_matrix_approx_eq squeeze(exp_obs[:startdate][1, :, :, :],    1) obs
# @test_matrix_approx_eq squeeze(exp_pseudo[:startdate][1, :, :, :], 1) pseudo

# With shockdec_startdate null
m <= Setting(:shockdec_startdate, Nullable{Date}())
states, obs, pseudo = shock_decompositions(m, system, histshocks)

# @test_matrix_approx_eq squeeze(exp_states[:no_startdate][1, :, :, :], 1) states
# @test_matrix_approx_eq squeeze(exp_obs[:no_startdate][1, :, :, :],    1) obs
# @test_matrix_approx_eq squeeze(exp_pseudo[:no_startdate][1, :, :, :], 1) pseudo


nothing
