#=
    This script runs the full estimation step with a preoptimized mode and a precomputed hessian, with precomputed
    random numbers.
    
    08/27/2015 ELM
=#
using DSGE

tic()
# Run without reoptimizing anything
m = Model990()
m.savepath = normpath(joinpath(dirname(@__FILE__),"../../save/m990-no_reoptimize_no_recalc_hessian/"))
createSaveDirectories(m.savepath)

# Since we want to run estimate with precalculated random vectors,
# reset the metropolis-hastings specifications for the whole simulation
m.num_mh_simulations_test= m.num_mh_simulations
m.num_mh_blocks_test = m.num_mh_blocks
m.num_mh_burn_test=m.num_mh_burn

estimate(m,verbose=true,testing=true,using_matlab_sigscale=true)
computeMoments(m)
time_elapsed = toq()

println("Time without recalculating Hessian: $time_elapsed")

