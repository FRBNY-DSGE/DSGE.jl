using DSGE

# construct a model object
m = Model990()
save = saveroot(m);
specify_mode!(m, "$save/output_data/m990/ss2/estimate/raw/paramsmode_vint=150102_iris.h5")
m <= Setting(:reoptimize, true)
#specify_hessian(m, "$save/output_data/m990/ss2/estimate/raw/hessian_vint=151127.h5")
#m <= Setting(:calculate_hessian, true)

# estimate as of 2015-Q3 using the default data vintage from 2015 Nov 27
m <= Setting(:data_vintage, "150102")
#m <= Setting(:data_vintage, "151127")
m <= Setting(:date_mainsample_end, quartertodate("2014-Q4"))

# reoptimize parameter vector, compute Hessian at mode, and full posterior
# parameter sampling
estimate(m;verbose=:low)

# produce LaTeX tables of parameter moments
compute_moments(m)
