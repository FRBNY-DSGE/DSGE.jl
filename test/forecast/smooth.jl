using DSGE, DataFrames, JLD, DistributedArrays
include("../util.jl")

path = dirname(@__FILE__())

# Set up arguments
custom_settings = Dict{Symbol, Setting}(
    :n_anticipated_shocks => Setting(:n_anticipated_shocks, 6),
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")),
    :use_parallel_workers => Setting(:use_parallel_workers, true))
m = Model990(custom_settings = custom_settings, testing = true)

df, system, kal = jldopen("$path/../reference/forecast_args.jld","r") do file
    read(file, "df"), read(file, "system"), read(file, "kal")
end

# Read expected output
exp_states, exp_shocks, exp_pseudo = jldopen("$path/../reference/smooth_out.jld", "r") do file
    read(file, "exp_states"),
    read(file, "exp_shocks"),
    read(file, "exp_pseudo")
end

# Call smoother and test
for smoother in [:durbin_koopman, :kalman]
    m <= Setting(:forecast_smoother, smoother)

    states, shocks, pseudo = smooth(m, df, system, kal)

    @test_matrix_approx_eq squeeze(exp_states[smoother, :z0][1, :, :], 1) states
    @test_matrix_approx_eq squeeze(exp_shocks[smoother, :z0][1, :, :], 1) shocks
    @test_matrix_approx_eq squeeze(exp_pseudo[smoother, :z0][1, :, :], 1) pseudo
end

nothing