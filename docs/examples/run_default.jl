using DSGE

# construct a model object
m = Model990()

# estimate as of 2015-Q3 using the default data vintage from 2015 Nov 27
m <= Setting(:data_vintage, "151127")
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))

# reoptimize parameter vector, compute Hessian at mode, and full posterior
# parameter sampling
estimate(m)

# produce LaTeX tables of parameter moments
compute_moments(m)

# forecast full distribution using 10 processes
my_procs = addprocs(10)
@everywhere using DSGE
m <= Setting(:use_parallel_workers, true)

# compute smoothed historical states, forecasted states and observables
output_vars = [:histstates, :forecaststates, :forecastobs]
forecast_one(m, :full, :none, output_vars)

# compute means and bands
compute_meansbands(m, :full, :none, output_vars)

rmprocs(my_procs)
