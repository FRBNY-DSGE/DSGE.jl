using Base.Test, DSGE

# Set up arguments
custom_settings = Dict{Symbol, Setting}(
    :n_anticipated_shocks => Setting(:n_anticipated_shocks, 6),
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2016-Q1")),
    :use_parallel_workers => Setting(:use_parallel_workers, true),
    :forecast_kill_shocks => Setting(:forecast_kill_shocks, true),
    :forecast_enforce_zlb => Setting(:forecast_enforce_zlb, false))
m = Model990(custom_settings = custom_settings, testing = true)

nstates = n_states_augmented(m)
nshocks = n_shocks_exogenous(m)
npseudo = n_pseudoobservables(m)
horizon = forecast_horizons(m)

system = compute_system(m)
z0 = zeros(nstates)

### 1. Supplying shocks

shocks = rand(nshocks, horizon)
states, obs, pseudo, shocks = compute_forecast(m, system, z0; shocks = shocks)
states, obs, pseudo, shocks = compute_forecast(system, z0, shocks)

# Enforce ZLB
ind_r = m.observables[:obs_nominalrate]
ind_r_sh = m.exogenous_shocks[:rm_sh]
zlb_value = forecast_zlb_value(m)
shocks = zeros(nshocks, horizon)
shocks[ind_r_sh, :] = -10.

states, obs, pseudo, shocks = compute_forecast(m, system, z0; shocks = shocks)
@assert all(x -> x < zlb_value, obs[ind_r, :])

m <= Setting(:forecast_enforce_zlb, true)
states, obs, pseudo, shocks = compute_forecast(m, system, z0; shocks = shocks)
@assert all(x -> abs(x - zlb_value) < 0.01, obs[ind_r, :])
@assert all(x -> x != -10.,                 shocks[ind_r_sh, :])
m <= Setting(:forecast_enforce_zlb, false)


### 2. Shocks not supplied => draw all zeros

states, obs, pseudo, shocks = compute_forecast(m, system, z0)
@assert all(x -> x == 0, shocks)


### 3. Drawing shocks

m <= Setting(:forecast_kill_shocks, false)

# Normally distributed shocks
states, obs, pseudo, shocks = compute_forecast(m, system, z0)
m <= Setting(:forecast_kill_shocks, true)

# t-distributed shocks
m <= Setting(:forecast_tdist_shocks, true)
states, obs, pseudo, shocks = compute_forecast(m, system, z0)
m <= Setting(:forecast_tdist_shocks, false)

m <= Setting(:forecast_kill_shocks, true)


nothing
