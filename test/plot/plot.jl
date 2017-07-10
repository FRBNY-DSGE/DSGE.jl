using DSGE, HDF5, JLD

path = dirname(@__FILE__)

# Initialize model object
m = AnSchorfheide(testing = true)
m <= Setting(:saveroot, tempdir())
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m <= Setting(:date_conditional_end, quartertodate("2015-Q4"))
m <= Setting(:use_population_forecast, true)
m <= Setting(:forecast_horizons, 8)

# Run full-distribution forecast
estroot = normpath(joinpath(dirname(@__FILE__), "..", "reference"))
overrides = forecast_input_file_overrides(m)
overrides[:full] = joinpath(estroot, "metropolis_hastings.h5")
output_vars = add_requisite_output_vars([:histobs, :forecastobs, :irfobs])

@everywhere using DSGE
m <= Setting(:forecast_block_size, 5)
@time forecast_one(m, :full, :none, output_vars, verbose = :none)
@time means_bands_all(m, :full, :none, output_vars; verbose = :none)

files = get_meansbands_output_files(m, :full, :none, output_vars)

# Plot history and forecast
hist_mb  = read_mb(files[:histobs])
fcast_mb = read_mb(files[:bddforecastobs])
output_file = joinpath(saveroot(m), "forecast__obs_nominalrate.pdf")
plot_history_and_forecast(:obs_nominalrate, hist_mb, fcast_mb;
                          output_file = output_file,
                          start_date = Nullable(DSGE.quartertodate("2007-Q1")))

# Plot forecast comparison
output_file = joinpath(saveroot(m), "forecastcomp__obs_nominalrate.pdf")
hist_mb2  = read_mb(files[:histobs])
fcast_mb2 = read_mb(files[:bddforecastobs])
plot_forecast_comparison(:obs_nominalrate, hist_mb, fcast_mb, hist_mb2, fcast_mb2;
                         output_file = output_file,
                         start_date = Nullable(DSGE.quartertodate("2007-Q1")))

# Hair plot
df = load_data(m)
output_file = joinpath(saveroot(m), "hairplot__obs_nominalrate.pdf")
hair_plot(:obs_nominalrate, df, [hist_mb], [fcast_mb];
          output_file = output_file)

# Plot IRF
irf_mb = read_mb(files[:irfobs])
output_file = joinpath(saveroot(m), "irf__rm_sh.pdf")
plot_irfs(:rm_sh, collect(keys(m.observables)), irf_mb;
          output_file = output_file)