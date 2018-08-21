using DSGE
using BenchmarkTools, JLD

path = dirname(@__FILE__)
filepath = "$path/../../reference/BondLabor/steady_state.jld"

# Benchmark current function
m = BondLabor()
steadystate!(m)
trial = @benchmark steadystate!($m) gcsample = true

# # Optionally over-write the existing reference trial
# write_ref_trial(trial, "steady_state", filepath = filepath)

# Compute time differential
print_all_benchmarks(trial, filepath, "steady_state")
