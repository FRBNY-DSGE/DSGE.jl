using DSGE
using BenchmarkTools, JLD

# Benchmark current function
m = KrusellSmith()
klein(m)
trial = @benchmark klein($m) gcsample = true

# Compute time differential
print_all_benchmarks(trial, "ref/klein_solve.jld")

# Optionally over-write the existing reference trial
# write_ref_trial(trial, "klein_solve")
