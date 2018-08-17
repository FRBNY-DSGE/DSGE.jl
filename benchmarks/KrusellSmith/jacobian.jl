using DSGE
using BenchmarkTools, JLD
import DSGE: jacobian

# Benchmark current function
m = KrusellSmith()
jacobian(m)
trial = @benchmark jacobian($m) gcsample = true

# Compute time differential
print_all_benchmarks(trial, "../reference/jacobian.jld", "jacobian")

# Optionally over-write the existing reference trial
# write_ref_trial(trial, "jacobian")
