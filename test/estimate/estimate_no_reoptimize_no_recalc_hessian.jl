#=
    This script runs the full estimation step with a preoptimized mode and a precomputed hessian,
    without precomputed random numbers.
    
    09/16/2015 ELM
=#
using DSGE

tic()

# Run without reoptimizing anything
m = Model990()
new_savepath = "/data/dsge_data_dir/dsgejl/estimate/save_noreop_norecalc/"

# Create directories
create_save_directories(m, new_savepath, reset_inpath=false) 

estimate(m,verbose=true,testing=false,using_matlab_sigscale=true)
time_elapsed = toq()

tic()
compute_moments(m)
time_elapsed_tables = toq()

println("Estimate, no reoptimize, no recalc Hessian, without using randvecs and randvals.")
println("Time without recalculating Hessian: $time_elapsed")
println("Time to compute tables = $time_elapsed_tables")
