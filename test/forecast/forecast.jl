using DSGE, JLD
include("../util.jl")

path = dirname(@__FILE__)

# Set up arguments
m = AnSchorfheide(testing = true)
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m <= Setting(:forecast_horizons, 1)

system, kal = jldopen("$path/../reference/forecast_args.jld","r") do file
    read(file, "system"), read(file, "kal")
end
z0 = zeros(n_states_augmented(m))

# Read expected output
exp_states, exp_obs, exp_pseudo, exp_shocks =
    jldopen("$path/../reference/forecast_out.jld", "r") do file
        read(file, "exp_states"),
        read(file, "exp_obs"),
        read(file, "exp_pseudo"),
        read(file, "exp_shocks")
    end

# Without shocks
states, obs, pseudo, shocks = forecast(m, system, z0; draw_shocks = false)

@test_matrix_approx_eq exp_states states
@test_matrix_approx_eq exp_obs    obs
@test_matrix_approx_eq exp_pseudo pseudo
@test_matrix_approx_eq exp_shocks shocks

# Supplying shocks
states, obs, pseudo, shocks = forecast(m, system, z0; shocks = shocks)

@test_matrix_approx_eq exp_states states
@test_matrix_approx_eq exp_obs    obs
@test_matrix_approx_eq exp_pseudo pseudo
@test_matrix_approx_eq exp_shocks shocks

# Draw normally distributed shocks
states, obs, pseudo, shocks = forecast(m, system, z0; draw_shocks = true)

# Draw t-distributed shocks
m <= Setting(:forecast_tdist_shocks, true)
states, obs, pseudo, shocks = forecast(m, system, z0; draw_shocks = true)
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


nothing