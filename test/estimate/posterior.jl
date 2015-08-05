using MATLAB
path = dirname(@__FILE__)

mf = MatFile("$path/posterior.mat")
YY = get_variable(mf, "YY")
lh_expected = get_variable(mf, "lnpy")
post_expected = get_variable(mf, "obj")
close(mf)

model = Model990()

lh = likelihood(model, YY)
@test_approx_eq lh_expected lh

post = posterior(model, YY)
@test_approx_eq post_expected post
