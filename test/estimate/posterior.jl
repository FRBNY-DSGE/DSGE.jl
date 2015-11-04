using DSGE
using HDF5, Base.Test

path = dirname(@__FILE__)

h5 = h5open("$path/../reference/posterior.h5")
YY = read(h5, "YY")
lh_expected = read(h5, "lnpy")
post_expected = read(h5, "obj")
close(h5)

model = Model990()

lh = likelihood(model, YY)
@test_approx_eq lh_expected lh

post = posterior(model, YY)
@test_approx_eq post_expected post

x = map(α->α.value, model.parameters)
post_at_start = posterior!(model, x, YY)
@test_approx_eq post_expected post_at_start

# Ensure if we are not evaluating at start vector, then we do not get the reference
# posterior
x = x .+ 0.01
post_not_at_start = posterior!(model, x, YY)
ϵ = 1.0
@test abs(post_at_start - post_not_at_start) > ϵ
