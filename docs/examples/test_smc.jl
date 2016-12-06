using DSGE

# Initialize a new instance of an An Schorfheide model 
m = AnSchorfheide()

# SMC can take advantage of parallel computing resources, should
# they be available
m <= Setting(:use_parallel_workers, false)

# SMC related settings
m <= Setting(:n_particles, 1000)
m <= Setting(:n_Φ, 200) # number of stages in SMC
m <= Setting(:λ, 2.0) # the curvature of the tempering schedule

# No need to find the mode for SMC
# If set to true, the mode will be calculated during estimation
m <= Setting(:reoptimize, false)
 
println("estimation beginning")
# If verbose=:low instead, the mean and sd of the model's
# parameters will not be displayed at each stage of SMC
estimate(m,verbose=:high,method=:SMC) 
println("estimation finished")

# Compute and save a table of moments
compute_moments(m)
println("tex table of moments saved in $(tablespath(m,"estimate"))")

