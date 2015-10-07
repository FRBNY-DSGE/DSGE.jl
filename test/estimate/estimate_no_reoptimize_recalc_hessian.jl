#=
    This script runs the estimate step where Julia reads in a pre-optimized mode but recalculates the Hessian. 
    
    09/17/2015 ELM
=#
using DSGE



# Run without reoptimizing anything
m = Model990()
new_savepath = "/data/dsge_data_dir/dsgejl/estimate/save_noreop_recalc/"
create_save_directories(m, new_savepath, reset_inpath=false)

m.num_mh_simulations_test= m.num_mh_simulations
m.num_mh_blocks_test = m.num_mh_blocks
m.num_mh_burn_test=m.num_mh_burn
m.recalculate_hessian = true

tic()
estimate(m,verbose=:low,testing=false)
time_elapsed = toq()

tic()
compute_moments(m)
time_elapsed_tables = toq()

println("Julia recalculates the hessian but reads in a pre-optimized mode.")
println("Time with recalculating Hessian: $time_elapsed")
println("Time for computing tables: $time_elapsed_tables")
