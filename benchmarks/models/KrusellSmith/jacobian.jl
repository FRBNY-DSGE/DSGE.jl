using DSGE
using BenchmarkTools, JLD
import DSGE: jacobian

path = dirname(@__FILE__)
filepath = "$path/../../reference/KrusellSmith/jacobian.jld"

# Benchmark current function
m = KrusellSmith()
jacobian(m)
trial = @benchmark jacobian($m) gcsample = true

# Optionally over-write the existing reference trial
# write_ref_trial(trial, "jacobian", filepath = filepath)

# Compute time differential
print_all_benchmarks(trial, filepath, "jacobian")
