using DSGE, ModelConstructors, Distributed
using Nullables, DataFrames, OrderedCollections, Dates

#############################
# To use:
# Just run in the Julia REPL
# include("run_default.jl")
# Note that the estimation
# step will take 2-3 hours.
#############################

##############
# Model Setup
##############
# Instantiate the FRBNY DSGE model object
m = Model1002("ss10")

# estimate as of 2015-Q3 using the default data vintage from 2015 Nov 27
m <= Setting(:data_vintage, "181115")
m <= Setting(:date_forecast_start, quartertodate("2018-Q4"))

# The following settings ensure that this script runs in
# a short amount of time. To properly estimate and
# forecast, we recommend either using the default settings
# (i.e. comment out the settings below) or
# changing the settings yourself.
m <= Setting(:n_mh_simulations, 100) # Do 100 MH steps during estimation
m <= Setting(:use_population_forecast, false) # Population forecast not available as data to turn off
m <= Setting(:forecast_block_size, 5) # adjust block size to run on small number of estimations

#############
# Estimation
#############
# Reoptimize parameter vector, compute Hessian at mode, and full posterior
# parameter sampling.
#
# Note some columns will have missing data because not all our data
# is publicly available. By default, `load_data` will error
# if any of the columns in the loaded DataFrame is empty,
# but we turn this feature off by setting the keyword `check_empty_columns = false`.
# Warnings will be still thrown after calling load_data indicating which columns
# are empty. However, estimate will still run when data is missing.
@time begin
    df = load_data(m, try_disk = false, check_empty_columns = false, summary_statistics = :none)
    data = df_to_matrix(m, df)
    estimate(m, data)
end

# produce LaTeX tables of parameter moments
moment_tables(m)

############
# Forecasts
############
# (1) forecast_one produces forecasts of the variables specified in output_vars
# and stores the forecasts in model units (log deviations from steady state).
# (2) compute_meansbands stores the forecast in a MeansBands type, defined in DSGE.jl,
# which transforms the forecast from model units to "human readable" units, i.e.
# the GDP forecast is transformed from annualized quarterly log per-capita growth rates
# to annualized quarterly aggregate percent change.

# The variables we want to forecast. In this case, all of the model observables
output_vars = [:histobs, :forecastobs]

# Modal forecast (point forecast)
forecast_one(m, :mode, :none, output_vars; check_empty_columns = false)
compute_meansbands(m, :mode, :none, output_vars; check_empty_columns = false)

# Full-distribution forecast (point forecast (mean) and uncertainty bands)
# Optionally add 10 processes to run the forecast in parallel (uncomment the 3 lines below).
# Alternatively, you can load the ClusterManagers package and add processes
# using one of the schedulers such as SGE or Slurm.
# addprocs(10)
# @everywhere using DSGE
# m <= Setting(:use_parallel_workers, true)

forecast_one(m, :full, :none, output_vars; check_empty_columns = false)
compute_meansbands(m, :full, :none, output_vars; check_empty_columns = false)

# Comment out the line below if you did not run the forecast in parallel.
# rmprocs(procs())
