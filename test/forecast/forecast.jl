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

system, kal = jldopen("$path/../reference/forecast_args.jld","r") do file
    read(file, "system"), read(file, "kal")
end
z0 = zeros(n_states_augmented(m))

# Read expected output
exp_states, exp_obs, exp_pseudo, exp_shocks = jldopen("$path/../reference/forecast_out.jld", "r") do file
    read(file, "exp_states"), read(file, "exp_obs"), read(file, "exp_pseudo"), read(file, "exp_shocks")
end

# Run forecast without supplying shocks
states, obs, pseudo, shocks = forecast(m, system, z0)

@test_matrix_approx_eq squeeze(exp_states[1, :, :], 1) states
@test_matrix_approx_eq squeeze(exp_obs[1, :, :],    1) obs
@test_matrix_approx_eq squeeze(exp_pseudo[1, :, :], 1) pseudo
@test_matrix_approx_eq squeeze(exp_shocks[1, :, :], 1) shocks

# Run forecast, supplying shocks
states, obs, pseudo, shocks = forecast(m, system, z0; shocks = shocks)

@test_matrix_approx_eq squeeze(exp_states[1, :, :], 1) states
@test_matrix_approx_eq squeeze(exp_obs[1, :, :],    1) obs
@test_matrix_approx_eq squeeze(exp_pseudo[1, :, :], 1) pseudo
@test_matrix_approx_eq squeeze(exp_shocks[1, :, :], 1) shocks

# Normally distributed shocks
m <= Setting(:forecast_kill_shocks, false)
states, obs, pseudo, shocks = forecast(m, system, z0)

# t-distributed shocks
m <= Setting(:forecast_tdist_shocks, true)
states, obs, pseudo, shocks = forecast(m, system, z0)
m <= Setting(:forecast_tdist_shocks, false)

# Enforce ZLB
ind_r = m.observables[:obs_nominalrate]
ind_r_sh = m.exogenous_shocks[:rm_sh]
zlb_value = forecast_zlb_value(m)
shocks = zeros(n_shocks_exogenous(m), forecast_horizons(m))
shocks[ind_r_sh, :] = -10.

states, obs, pseudo, shocks = forecast(m, system, z0; shocks = shocks)
@assert all(x -> x < zlb_value, obs[ind_r, :])

states, obs, pseudo, shocks = forecast(m, system, z0; shocks = shocks, enforce_zlb = true)
@assert all(x -> abs(x - zlb_value) < 0.01, obs[ind_r, :])
@assert all(x -> x != -10.,                 shocks[ind_r_sh, :])

# Draw z0
m <= Setting(:forecast_draw_z0, true)
states, obs, pseudo, shocks = forecast(m, system, kal)


nothing
