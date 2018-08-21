using DSGE
using BenchmarkTools, JLD

path = dirname(@__FILE__)
filepath = "$path/../../reference/BondLabor/klein.jld"

# Benchmark current function
m = BondLabor()
steadystate!(m)
trial = @benchmark klein($m) gcsample = true

# # Optionally over-write the existing reference trial
# write_ref_trial(trial, "klein", filepath = filepath)

# Compute time differential
print_all_benchmarks(trial, filepath, "klein")
