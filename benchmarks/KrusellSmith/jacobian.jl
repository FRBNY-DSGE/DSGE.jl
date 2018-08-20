using DSGE
using BenchmarkTools, JLD
import DSGE: jacobian

path = dirname(@__FILE__)

# Benchmark current function
m = KrusellSmith()
jacobian(m)
trial = @benchmark jacobian($m) gcsample = true

# Compute time differential
print_all_benchmarks(trial, "$path/../reference/jacobian.jld", "jacobian")

# Optionally over-write the existing reference trial
# write_ref_trial(trial, "jacobian")
