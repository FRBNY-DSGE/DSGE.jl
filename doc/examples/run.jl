using DSGE

# construct a model object
m = Model990()
save = saveroot(m);
specify_mode!(m, "$save/output_data/m990/ss2/estimate/raw/paramsmode_vint=150102_iris.h5")
#m <= Setting(:reoptimize, true)
#m <= Setting(:n_hessian_test_params, 3)

specify_hessian(m, "$save/output_data/m990/ss2/estimate/raw/hessian_vint=150102_iris.h5")
#m <= Setting(:calculate_hessian, true)

# estimate as of 2015-Q3 using the default data vintage from 2015 Nov 27
m <= Setting(:data_vintage, "150102")
#m <= Setting(:data_vintage, "151127")
m <= Setting(:date_mainsample_end, quartertodate("2014-Q4"))

# reoptimize parameter vector, compute Hessian at mode, and full posterior
# parameter sampling
#file = h5open("$save/output_data/m990/ss2/estimate/raw/propcov.h5","r")
#propdist_cov  = read(file, "propcov")
#estimate(m;verbose=:high, proposal_covariance = propdist_cov)
estimate(m;verbose=:high)

# produce LaTeX tables of parameter moments
compute_moments(m)
