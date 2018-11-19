using DSGE
using Nullables, DataFrames, OrderedCollections, Dates

##############
# Model Setup
##############
# Instantiate the FRBNY DSGE model object
m = Model1002()

# estimate as of 2015-Q3 using the default data vintage from 2015 Nov 27
m <= Setting(:data_vintage, "181115")
m <= Setting(:date_forecast_start, quartertodate("2018-Q4"))

#############
# Estimation
#############
# reoptimize parameter vector, compute Hessian at mode, and full posterior
# parameter sampling
estimate(m)

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

m <= Setting(:use_population_forecast, false)

# Modal forecast (point forecast)
forecast_one(m, :mode, :none, output_vars)
compute_meansbands(m, :mode, :none, output_vars)

# Full-distribution forecast (point forecast (mean) and uncertainty bands)
# Optionally add 10 processes to run the forecast in parallel (uncomment the 3 lines below).
# Alternatively, you can load the ClusterManagers package and add processes
# using one of the schedulers such as SGE or Slurm.
addprocs(10)
@everywhere using DSGE
m <= Setting(:use_parallel_workers, true)

forecast_one(m, :full, :none, output_vars)
compute_meansbands(m, :full, :none, output_vars)

rmprocs(procs())
