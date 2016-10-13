using DSGE
# using Distributions

# Set up arguments
custom_settings = Dict{Symbol, Setting}(
    :n_anticipated_shocks => Setting(:n_anticipated_shocks, 6),
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2016-Q1")),
    :use_parallel_workers => Setting(:use_parallel_workers, true),
    :forecast_kill_shocks => Setting(:forecast_kill_shocks, true))
m = Model990(custom_settings = custom_settings, testing = true)

nstates = n_states_augmented(m)
nshocks = n_shocks_exogenous(m)
npseudo = n_pseudoobservables(m)
horizon = forecast_horizons(m)

system = compute_system(m)
z0 = zeros(nstates)

# Test invocations supplying shocks
shocks = rand(nshocks, horizon)
states, obs, pseudo, shocks = compute_forecast(m, system, z0; shocks = shocks)
states, obs, pseudo, shocks = compute_forecast(system, z0, shocks)

# Test invocations not supplying shocks
states, obs, pseudo, shocks = compute_forecast(m, system, z0)

m <= Setting(:forecast_tdist_shocks, true)
states, obs, pseudo, shocks = compute_forecast(m, system, z0)

m <= Setting(:forecast_kill_shocks, true)
states, obs, pseudo, shocks = compute_forecast(m, system, z0)

nothing
