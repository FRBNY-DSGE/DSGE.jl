using DSGE
using HDF5, Base.Test

path = dirname(@__FILE__)

m = AnSchorfheide()
m.testing = true

file = "$path/../reference/posterior.h5"
data = h5read(file, "data")
lh_expected = h5read(file, "likelihood")
post_expected = h5read(file, "posterior")

lh = likelihood(m, data)
@test_approx_eq lh_expected lh

post = posterior(m, data)
@test_approx_eq post_expected post

x = map(α->α.value, m.parameters)
post_at_start = posterior!(m, x, data)
@test_approx_eq post_expected post_at_start

# Ensure if we are not evaluating at start vector, then we do not get the reference
# posterior
y = x .+ 0.01
post_not_at_start = posterior!(m, y, data)
ϵ = 1.0
@test abs(post_at_start - post_not_at_start) > ϵ