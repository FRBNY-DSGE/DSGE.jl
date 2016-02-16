using DSGE
using HDF5, Base.Test

path = dirname(@__FILE__)

h5 = h5open("$path/../reference/posterior.h5")
data = read(h5, "data")
lh_expected = read(h5, "lnpy")
post_expected = read(h5, "obj")
close(h5)

m = Model990(testing=true)

lh, _ = likelihood(m, data)
@test_approx_eq lh_expected lh

post = posterior(m, data)[:post]
@test_approx_eq post_expected post

x = map(α->α.value, m.parameters)
post_at_start = posterior!(m, x, data)[:post]
@test_approx_eq post_expected post_at_start

# Ensure if we are not evaluating at start vector, then we do not get the reference
# posterior
y = x .+ 0.01
post_not_at_start = posterior!(m, y, data)[:post]
ϵ = 1.0
@test abs(post_at_start - post_not_at_start) > ϵ
