using DSGE, HDF5, JLD

path = dirname(@__FILE__)

# Initialize model object
m = AnSchorfheide(testing = true)
m <= Setting(:saveroot, tempdir())
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m <= Setting(:date_conditional_end, quartertodate("2015-Q4"))
m <= Setting(:use_population_forecast, true)
m <= Setting(:forecast_horizons, 8)
m <= Setting(:impulse_response_horizons, 8)

# Run full-distribution forecast
estroot = normpath(joinpath(dirname(@__FILE__), "..", "reference"))
overrides = forecast_input_file_overrides(m)
overrides[:full] = joinpath(estroot, "metropolis_hastings.h5")
output_vars = add_requisite_output_vars([:histobs, :forecastobs, :shockdecobs, :irfobs])

@everywhere using DSGE
m <= Setting(:forecast_block_size, 5)
@time forecast_one(m, :full, :none, output_vars, verbose = :none)
@time means_bands_all(m, :full, :none, output_vars; verbose = :none)

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

# Hair plot
df = load_data(m)
output_file = joinpath(saveroot(m), "hairplot__obs_nominalrate.pdf")
hist_mb = read_mb(m, :full, :none, :histobs)
fcast_mb = read_mb(m, :full, :none, :bddforecastobs)
hair_plot(:obs_nominalrate, df, [hist_mb], [fcast_mb];
          output_file = output_file, verbose = :none)