using DSGE
using HDF5

# Initialize a new instance of an An Schorfheide model
m = AnSchorfheide()

# Data can be constructed automatically, i.e. calling estimate(m),
# but for ease of example we will load in data
data = h5read("../../test/reference/smc.h5", "data")

# SMC can take advantage of parallel computing resources, should
# they be available
m <= Setting(:use_parallel_workers, false)

# SMC related settings
m <= Setting(:n_particles, 800)
m <= Setting(:n_Φ, 50) # number of stages in SMC
m <= Setting(:λ, 2.0) # the curvature of the tempering schedule

# No need to find the mode for SMC
# If set to true, the mode will be calculated during estimation
m <= Setting(:reoptimize, false)
m <= Setting(:sampling_method, :SMC)

println("estimation beginning")
# If verbose=:low instead, the mean and sd of the model's
# parameters will not be displayed at each stage of SMC
estimate(m, data, verbose = :high)
println("estimation finished")
println("sample saved in $(rawpath(m, "estimate", "smcsave.h5")))")
println("particle cloud saved in $(rawpath(m, "estimate", "smc_cloud.jld")))")

# Compute and save a table of moments
moment_tables(m)
println("tex table of moments saved in $(tablespath(m, "estimate"))")

