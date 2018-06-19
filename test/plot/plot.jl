using DSGE, HDF5, JLD, Plots

path = dirname(@__FILE__)

# Initialize the plotting backend
gr()

# Initialize model object
m = AnSchorfheide(testing = true)
m <= Setting(:saveroot, tempdir())
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m <= Setting(:date_conditional_end, quartertodate("2015-Q4"))
m <= Setting(:use_population_forecast, true)
m <= Setting(:forecast_horizons, 12)
m <= Setting(:impulse_response_horizons, 8)

# Run full-distribution forecast
estroot = normpath(joinpath(dirname(@__FILE__), "..", "reference"))
overrides = forecast_input_file_overrides(m)
overrides[:mode] = joinpath(estroot, "optimize.h5")
overrides[:full] = joinpath(estroot, "metropolis_hastings.h5")
output_vars = add_requisite_output_vars([:histobs, :forecastobs, :shockdecobs, :irfobs])

@everywhere using DSGE
m <= Setting(:forecast_block_size, 5)
@time forecast_one(m, :full, :none, output_vars, verbose = :none)
@time compute_meansbands(m, :full, :none, output_vars; verbose = :none)

# Plot history and forecast
plot_history_and_forecast(m, :obs_nominalrate, :obs, :full, :none,
                          bdd_and_unbdd = true,
                          start_date = DSGE.quartertodate("2007-Q1"),
                          verbose = :none)

# Plot forecast comparison
plot_forecast_comparison(m, m, :obs_nominalrate, :obs, :full, :none,
                         bdd_and_unbdd = true,
                         start_date = DSGE.quartertodate("2007-Q1"),
                         verbose = :none)

# Plot shock decomposition
plot_shock_decomposition(m, :obs_nominalrate, :obs, :full, :none, verbose = :none)

# Plot IRF
plot_impulse_response(m, :rm_sh, collect(keys(m.observables)), :obs, :full, :none,
                      verbose = :none)
plot_impulse_response(m, :rm_sh, collect(keys(m.observables)), :obs, :full, :none,
                      flip = true, verbose = :none)

# Plot scenario, zeroing out measurement error
for para in [:e_y, :e_π, :e_R]
    m[para].value = 0
end
alt = Scenario(:altscen, "Test Alternative Scenario", [:obs_gdp, :obs_cpi], [:g_sh, :rm_sh], "REF")
forecast_scenario(m, alt, verbose = :none)
scenario_means_bands(m, alt, verbose = :none)
plot_scenario(m, :obs_nominalrate, :obs, alt, untrans = true, verbose = :none)
plot_scenario(m, :obs_nominalrate, :obs, alt, verbose = :none)
plot_scenario(m, :obs_nominalrate, :obs, alt, fourquarter = true, verbose = :none)
@test_throws ErrorException plot_scenario(m, :obs_nominalrate, :obs, alt,
                                          untrans = true, fourquarter = true)

# Hair plot
realized = load_data(m, verbose = :none)
hist_mb = read_mb(m, :full, :none, :histobs)
fcast_mb = read_mb(m, :full, :none, :bddforecastobs)
hair_plot(:obs_nominalrate, realized, [hist_mb], [fcast_mb];
          plotroot = saveroot(m), verbose = :none)

# Forecast decomposition
m_old = deepcopy(m)
m_old <= Setting(:date_forecast_start, quartertodate("2014-Q4"))
m_old <= Setting(:date_conditional_end, quartertodate("2014-Q4"))
df_new = load_data(m)
df_old = df_new[1:end-4, :]
@time decompose_forecast(m, m_old, df_new, df_old, :mode, :none, :none, [:obs]; verbose = :none)
@time decomposition_means(m, m_old, :mode, :none, :none, [:obs], verbose = :none)
for indshocks in [true, false]
    plot_forecast_decomposition(m, m_old, [:obs_nominalrate], :obs, :mode, :none, :none,
                                individual_shocks = indshocks, verbose = :none)
end