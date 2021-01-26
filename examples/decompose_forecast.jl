using DSGE, ModelConstructors, Distributed
using DataFrames, OrderedCollections, Dates, Plots
GR.inline("pdf")

####################################
# To use:
# Just run in the Julia REPL
# include("decompose_forecast.jl")
####################################

use_parallel = false # set to true if you want to use parallel workers
n_workers    = 10
make_plots   = true  # set to true if you want to produce forecast decomposition plots

############
# Estimate
############

# Initialize a new instance of an An Schorfheide model,
# a canonical New Keynesian DSGE model augmented for
# Bayesian estimation.
m = AnSchorfheide()

# See "test_smc.jl" for more information on setting up SMC
# We choose settings so the estimation finishes faster,
# so the particles will not approximate the actual posterior well.
m <= Setting(:data_vintage, "200101", true, "vint", "Data vintage")
m <= Setting(:date_forecast_start, DSGE.quartertodate("2020-Q1"))
m <= Setting(:test_decomp_fcast, "true", true, "test_decomp_fcast",
             "Tag for estimation output file name indicating this file is for testing a forecast decomposition")
m <= Setting(:use_parallel_workers, use_parallel)
m <= Setting(:n_particles, 800)
m <= Setting(:n_Φ, 50) # number of stages in SMC
m <= Setting(:λ, 2.0) # the curvature of the tempering schedule
m <= Setting(:adaptive_tempering_target_smc, false) # Setting a fixed tempering schedule
m <= Setting(:reoptimize, false)
m <= Setting(:sampling_method, :SMC)
forecast_string = "new" # string to distinguish the forecast for the new model

# Load data and create fake old data
df = load_data(m)
df_old = deepcopy(df[1:end - 4, :]) # subset data into a pseudo "old" and "new" data so old data is 1 year old

# Set up our "fake" old model
m_old = deepcopy(m)
m_old <= Setting(:data_vintage, "190101", true, "vint", "Data vintage") # fake data vintage, 1 year old
m_old <= Setting(:date_forecast_start, DSGE.quartertodate("2019-Q1"))
forecast_string_old = "old" # forecast string for the old forecast

# Estimate
if use_parallel
    myprocs = addprocs(n_workers)
    @everywhere using DSGE, OrderedCollections
end

estimate(m, df_to_matrix(m, df), verbose = :low, run_csminwel = false)
estimate(m_old, df_to_matrix(m_old, df_old), verbose = :low, run_csminwel = false)

if use_parallel
    rmprocs(myprocs)
end

#####################
# Decompose forecast
#####################
# We now produce a forecast decomposition, which decomposes the difference
# in the forecasts of the old model/data and the new model/data into the effect from
# data revisions, news (e.g. new data), and parameter re-estimation.
# To save on time, we only compute the modal forecast decomposition,
# and we will not use the DSGE's indivdual shocks for the decomposition.

# Run forecast
DSGE.decompose_forecast(m, m_old, df, df_old, :mode, :none, :none, [:obs, :pseudo],
                        check = true, forecast_string_new = forecast_string,
                        forecast_string_old = forecast_string_old)
DSGE.decomposition_means(m, m_old, :mode, :none, :none, [:obs, :pseudo],
                         verbose = :none, forecast_string_new = forecast_string,
                         forecast_string_old = forecast_string_old)

if make_plots
    # Set up plots
    kwargs = Dict{Symbol, Any}()
    kwargs[:tick_size]      = 1
    kwargs[:hist_label]     = ""
    kwargs[:forecast_label] = ""
    kwargs[:start_date]     = DSGE.iterate_quarters(date_mainsample_end(m), -4)
    kwargs[:end_date]       = DSGE.iterate_quarters(date_mainsample_end(m), 8)  # 2 years out

    # Make plots
    DSGE.plot_forecast_decomposition(m, m_old, collect(keys(m.observables)),
                                     :obs, :mode, :none, :none;
                                     individual_shocks = false,
                                     forecast_string_new = forecast_string,
                                     forecast_string_old = forecast_string_old,
                                     kwargs...)
    DSGE.plot_forecast_decomposition(m, m_old, collect(keys(m.pseudo_observables)),
                                     :pseudo, :mode, :none, :none;
                                     individual_shocks = false,
                                     forecast_string_new = forecast_string,
                                     forecast_string_old = forecast_string_old,
                                     kwargs...)
end

nothing
