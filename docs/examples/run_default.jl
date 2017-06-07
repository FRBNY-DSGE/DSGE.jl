using DSGE

# estimate as of 2015-Q3 using the default data vintage from 2015 Nov 27
custom_settings = Dict{Symbol, Setting}(
    :data_vintage        => Setting(:data_vintage, "151127"),
    :date_forecast_start => Setting(:date_forecast_start, quartertodate("2015-Q4")))

# construct a model object
m = Model990()

# reoptimize parameter vector, compute Hessian at mode, and full posterior
# parameter sampling
estimate(m)

# produce LaTeX tables of parameter moments
compute_moments(m)

# forecast using 10 processes
my_procs = addprocs(10)
@everywhere using DSGE

m <= Setting(:use_parallel_workers, true)
output_vars = get_output_vars(m, :simple)
forecast_all(m, [:none, :semi, :full], [:mode, :full], output_vars; procs = my_procs)
rmprocs(my_procs)
