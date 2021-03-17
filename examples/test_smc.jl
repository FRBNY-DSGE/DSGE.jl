using DSGE, ModelConstructors
using JLD2

# Initialize a new instance of an An Schorfheide model
m = AnSchorfheide()

# Data can be constructed automatically, i.e. calling estimate(m),
# but for ease of example we will load in data
@load "../test/reference/smc.jld2" data
data = Matrix{Float64}(data) # may be loaded as Matrix{Union{Float64, Missing}}, but shouldn't have missing data, so we recast the type here.

# SMC can take advantage of parallel computing resources, should
# they be available
m <= Setting(:use_parallel_workers, false)

# SMC related settings
m <= Setting(:n_particles, 800)
m <= Setting(:n_Φ, 50) # number of stages in SMC
m <= Setting(:λ, 2.0) # the curvature of the tempering schedule
m <= Setting(:adaptive_tempering_target_smc, false) # Setting a fixed tempering schedule

# No need to find the mode for SMC
# If set to true, the mode will be calculated during estimation
m <= Setting(:reoptimize, false)
m <= Setting(:sampling_method, :SMC)

println("estimation beginning")
# If verbose=:low instead, the mean and sd of the model's
# parameters will not be displayed at each stage of SMC
# By default, when using SMC, we run csminwel after
# the estimation completes to find the true posterior mode,
# but this can take several hours, even if SMC finishes relatively quickly.
estimate(m, data, verbose = :low, run_csminwel = false)
println("estimation finished")
println("sample saved in $(rawpath(m, "estimate", "smcsave.h5")))")
println("particle cloud saved in $(rawpath(m, "estimate", "smc_cloud.jld2")))")

# Compute and save a table of moments
moment_tables(m)
println("tex table of moments saved in $(tablespath(m, "estimate"))")
