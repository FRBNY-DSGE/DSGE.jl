using Base: Test

using HDF5
using DSGE
include("../util.jl")

m = Model990()
m.reoptimize = false
m.recalculate_hessian = false
m.num_mh_blocks_test = 1
m.num_mh_simulations_test = 100
m.num_mh_burn_test = 0
m.mh_thinning_step = 1

estimate(m, verbose=:none, testing=true, using_matlab_sigscale=true)

# Read in the parameter draws from estimate()
h5_fn = joinpath(outpath(m), "sim_save.h5")
h5 = h5open(h5_fn, "r")
jl_draws = read(h5, "parasim")
close(h5)

# Read in the reference parameter draws
reference_fn = joinpath(dirname(@__FILE__), "metropolis_hastings.h5")
h5_ref = h5open(reference_fn, "r")
ref_draws = read(h5_ref, "theta")
close(h5_ref)

# Test equal
@test test_matrix_eq(ref_draws, jl_draws, ε=1e-9)
