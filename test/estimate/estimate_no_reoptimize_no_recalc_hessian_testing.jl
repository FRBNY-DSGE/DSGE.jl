#=
    This script runs the full estimation step with a preoptimized mode and a precomputed hessian, with precomputed
    random numbers. It should correspond to the output from Main_norecalc_testing.m in ~/dsge/cleanCode990.
    
    09/16/2015 ELM
=#
using DSGE

# Run without reoptimizing anything
m = Model990()
#new_savepath = normpath(joinpath(dirname(@__FILE__),"../../save/m990-no_reoptimize_no_recalc_hessian/"))
new_savepath = "/data/dsge_data_dir/dsgejl/estimate/save_noreop_norecalc_testing/"

# Create directories
createSaveDirectories(m, new_savepath, reset_inpath=false) 

# Since we want to run estimate with precalculated random vectors,
# reset the metropolis-hastings specifications for the whole simulation
m.num_mh_simulations_test= m.num_mh_simulations
m.num_mh_blocks_test = m.num_mh_blocks
m.num_mh_burn_test=m.num_mh_burn

tic()
estimate(m,verbose=true,testing=true,using_matlab_sigscale=true)
time_elapsed = toq()

tic()
computeMoments(m)
time_elapsed_tables = toq()

println("No reoptimize, no recalc hessian, using precomputed random vectors and values.")
println("Time without recalculating Hessian: $time_elapsed")
println("Time to compute tables: $time_elapsed_tables")

