#=
Goal

To see if the output of a single step in the RWMH algorithm as it is implemented currently in
estimate.jl is equivalent to the output of a single step implemented in mutation_RWMH.jl

Procedure

Choose a set initial distribution, likelihood, and thus a posterior.
Set the tuning parameter equal to 1 (thus there is no tempering since tempering is not a 
characteristic of our current implementation of metropolis hastings.
Compute reference "mutations" from our estimate.jl metropolis hastings
Compute test "mutations" from mutation_RWMH.jl
Test that the proposed steps are equal i.e.
Test that the proposed/accepted or rejected step is the same for both algos
For best coverage test for both an accept and a reject
=#

using DSGE
using DSGEModels
using HDF5, Base.Test
using DataFrames

path = dirname(@__FILE__)
include("../util.jl")

# Set up models for testing
s = Schorf()
s.testing=true

s <= Setting(:saveroot, "$path/../../save")
s <= Setting(:dataroot, "$path/../../save/input_data")
s <= Setting(:data_vintage, "160718")
s <= Setting(:date_presample_start, quartertodate("1983-Q1"))
s <= Setting(:date_mainsample_end, quartertodate("2002-Q4"))

# Read in the data
df = load_data(s)

# Reading in reference stuff
file = h5open("$path/../reference/mutation_RWMH.h5","r")

acc_rvecs = file["Acc_rvecs"][:,:]
acc_rvals = file["Acc_rvals"][:]
acc_init_likelihoods = file["Acc_init_likelihoods"][:]
acc_init_para = file["Acc_init_para"][:,:]
acc_init_post = file["Acc_init_post"][:]
acc_prop_likelihoods = file["Acc_prop_likelihoods"][:]
acc_prop_para = file["Acc_prop_para"][:,:]
acc_prop_post = file["Acc_prop_post"][:]

rej_rvecs = file["Rej_rvecs"][:,:]
rej_rvals = file["Rej_rvals"][:]
rej_init_likelihoods = file["Rej_init_likelihoods"][:]
rej_init_para = file["Rej_init_para"][:,:]
rej_init_post = file["Rej_init_post"][:]
rej_prop_likelihoods = file["Rej_prop_likelihoods"][:]
rej_prop_para = file["Rej_prop_para"][:,:]
rej_prop_post = file["Rej_prop_post"][:]

tune_phi = file["tune_phi"][:]
tune_R = file["tune_R"][:,:]
i = 2 #choosing the first recursion of SMC for testing purposes

close(file)

type TestTune
    npara::Int64
    c::AbstractFloat
    phi::Array{Float64}
    R::AbstractArray 
end

tune = TestTune(13, 0.5, tune_phi, tune_R)

return

#Quick and dirty solution only testing one accept and one reject. 
#Think about a way to generalize this further to compare matrix output?
atest_para, atest_loglh, atest_post, atest_acpt = mutation_RWMH(acc_init_para[:,1], acc_init_likelihoods[1], acc_init_post[1], tune, i, df, s; rvec = acc_rvecs[1], rval = acc_rvals[1], px = acc_prop_para[:,1], lx = acc_prop_likelihoods[1], postx = acc_prop_post[1])

@test_approx_eq acc_prop_para[:,1] atest_para
@test_approx_eq acc_prop_likelihoods[1] atest_loglh
@test_approx_eq acc_prop_post[1] atest_post
@test_approx_eq atest_acpt 1

rtest_para, rtest_loglh, rtest_post, rtest_acpt = mutation_RWMH(rej_init_para[:,1], rej_init_likelihoods[1], rej_init_post[1], tune, i, df, s; rvec = rej_rvecs[1], rval = rej_rvals[1], px = rej_prop_para[:,1], lx = rej_prop_likelihoods[1], postx = rej_prop_post[1])

@test_approx_eq rej_init_para[:,1] rtest_para
@test_approx_eq rej_init_likelihoods[1] rtest_loglh
@test_approx_eq rej_init_post[1] rtest_post
@test_approx_eq rtest_acpt 0

nothing
