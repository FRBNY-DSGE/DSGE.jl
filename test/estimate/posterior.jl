using HDF5

path = dirname(@__FILE__)

h5 = h5open("$path/posterior.h5")
YY = read(h5, "YY")
lh_expected = read(h5, "lnpy")
post_expected = read(h5, "obj")
close(h5)

model = Model990()

lh = likelihood(model, YY)
@test_approx_eq lh_expected lh

post = posterior(model, YY)
@test_approx_eq post_expected post
