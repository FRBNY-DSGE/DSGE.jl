#=
    This is a skeleton of what will ultimately test the entire estimation routine. For now,
    be careful, as this will take >24 hours.
=#
using DSGE

m = Model990()
m.reoptimize=true
m.recalculate_hessian=true
estimate(m, verbose=true)
