#=
    This script runs the estimation step twice, both with a pre-optimized mode. The first
    runs with a pre-optimized Hessian, the second recalculates the Hessian before proceeding
    to Metropolis-Hastings.
    
    08/27/2015 ELM
=#
using DSGE


tic()
# Run without reoptimizing anything
m = Model990()
m.savepath = normpath(joinpath(dirname(@__FILE__),"../../save/m990-no_reoptimize_no_recalc_hessian/"))
createSaveDirectories(m.savepath)
estimate(m,verbose=true,testing=true,using_matlab_sigscale=true)
computeMoments(m)
time_elapsed = toq()

println("Time without recalculating Hessian: $time_elapsed")

