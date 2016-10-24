using DSGE

# construct a model object
m = Model990()

# estimate as of 2015-Q3 using the default data vintage from 2015 Nov 27
m <= Setting(:data_vintage, "151127")
m <= Setting(:date_mainsample_end, quartertodate("2015-Q3"))

# reoptimize parameter vector, compute Hessian at mode, and full posterior
# parameter sampling
estimate(m)

# produce LaTeX tables of parameter moments
compute_moments(m)
