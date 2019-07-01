using DSGE
using BenchmarkTools, JLD

path = dirname(@__FILE__)
filepath = "$path/../../reference/BondLabor/jacobian.jld"

# Benchmark current function
m = BondLabor()
steadystate!(m)
trial = @benchmark DSGE.jacobian($m) gcsample = true

# # Optionally over-write the existing reference trial
# write_ref_trial(trial, "jacobian", filepath = filepath)

# Compute time differential
print_all_benchmarks(trial, filepath, "jacobian")
