using ClusterManagers, DSGE, ModelConstructors, FileIO, Plots, Distributed
GR.inline("pdf")
fn = dirname(@__FILE__)

# This script auto-generates a packet of
# plots and figures to analyze
# estimation, forecast, and impulse response results.
# This script also serves as an example of how we recommend
# structuring specification or "spec" files for running forecasts, i.e.
# this file has a recommended workflow for the
# entirety of a forecast.

# A note also on the options. Generally, the only memory-intensive
# operation is the computation of the shock decompositions,
# which involves handling many large matrices. If memory
# needs to be pre-allocated, then 1GB is often sufficient
# if you do not run the shock decompositions portion of the code.

# What do you want to do?
run_full_forecast = true  # Run full distribution forecast
do_histforecast   = true  # Compute history and forecast of observables and pseudo-observables
do_shockdecs      = true  # Compute shock decompositions
do_irfs           = true  # Compute impulse response functions
make_packet       = true  # Create packet w/info on estimation, forecasts, and IRFs
make_plots        = true  # Make plots
add_workers       = false # Run in parallel
n_workers         = 10

# Initialize model objects and desired settings
m = Model1002("ss10")
m <= Setting(:sampling_method, :MH)
usual_model_settings!(m, "181115"; cdvt = "181115", fcast_date = quartertodate("2018-Q4"))
m <= Setting(:use_population_forecast, false) # We do not make this data available
m <= Setting(:saveroot, "$(fn)/../save/") # set save root
m <= Setting(:dataroot, "$(fn)/../save/input_data/") # set data root
m <= Setting(:date_forecast_end, quartertodate("2030-Q1")) # when to stop forecast
m <= Setting(:forecast_horizons,                           # number of quarters in forecast horizon; must match
             DSGE.subtract_quarters(date_forecast_end(m), date_forecast_start(m)) + 1) # date_forecast_end to run successfully
forecast_string = "" # to identify different forecasts from each other
m <= Setting(:forecast_block_size, 40) # this depends on the number of samples from estimations
# and thinning during forecast (determined by
# the setting :forecast_jstep)
if get_setting(m, :sampling_method) == :SMC
    # If SMC is the estimation method, then we need t adjust some default settings
    m <= Setting(:forecast_jstep, 1)
end

# Load data to create data set. There will be empty columns
# since not all the data used for Model1002 is available from FRED.
# This requires setting the check_empty_columns keyword to false.
df = load_data(m; check_empty_columns = false)

# Use overrides to give the correct file path for saved estimation output
overrides = forecast_input_file_overrides(m)
overrides[:full] = "$(fn)/../test/reference/mhsave_vint=181115.h5"

# Full-distribution forecast
if run_full_forecast

    output_vars = Vector{Symbol}(undef,0)
    if do_histforecast
        # Write data to create historical and forecast output
        output_vars = vcat(output_vars, [:histpseudo, :histobs, :histstdshocks,
                                         :hist4qpseudo, :hist4qobs, :histutpseudo,
                                         :forecastpseudo, :forecastobs, :forecastutpseudo,
                                         :forecast4qpseudo, :forecast4qobs, :forecaststdshocks])
    end

    if do_shockdecs
        # Shock decompositions of forecasts
        output_vars = vcat(output_vars, [:dettrendobs, :dettrendpseudo, :trendobs,
                                         :trendpseudo, :shockdecpseudo, :shockdecobs])
    end

    if do_irfs
        # Impulse response to all shocks, including for endogenous states
        output_vars = vcat(output_vars, [:irfobs, :irfpseudo]) # :irfstates not calculated by default b/c takes a long time
    end

    if add_workers
        my_procs = addprocs(n_workers)
        @everywhere using DSGE, OrderedCollections
    end

    usual_model_forecast(m, :full, :none, output_vars,
                         est_override = overrides[:full],
                         forecast_string = forecast_string,
                         density_bands = [.5, .6, .68, .7, .8, .9],
                         check_empty_columns = false)
end

# Make plots and packet
output_vars = [:forecastobs, :forecastpseudo]
if do_shockdecs
    output_vars = vcat(output_vars, [:shockdecobs, :shockdecpseudo])
end

if do_irfs
    # we do not plot irfs for endogenous states b/c too many!
    output_vars = vcat(output_vars, [:irfobs, :irfpseudo])
    sections = [:estimation, :forecast, :irf]
else
    sections = [:estimation, :forecast]
end

if make_plots
    plot_standard_model_packet(m, :full, :none, output_vars,
                               forecast_string = forecast_string,
                               sections = sections)
end

if make_packet
    # Choose forecast centered or the standard packet.
    # Note that these functions write to the same tex file,
    # so you can use only one for a given data vintage!
    # write_forecast_centric_packet(m, :full, :none, output_vars,
    #                               sections = sections, forecast_string = forecast_string)
    write_standard_model_packet(m, :full, :none, output_vars,
                                sections = sections, forecast_string = forecast_string)
    moment_tables(m)
end

if add_workers
    rmprocs(my_procs)
end

nothing
