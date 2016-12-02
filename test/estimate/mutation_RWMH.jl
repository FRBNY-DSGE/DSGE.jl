
using DSGE
using HDF5, Base.Test
using DataFrames

path = dirname(@__FILE__)
include("../util.jl")

# Set up models for testing
s = AnSchorfheide()
s.testing=true

s <= Setting(:saveroot, "$path/../../save")
s <= Setting(:dataroot, "$path/../../save/input_data")
s <= Setting(:data_vintage, "160706")
s <= Setting(:date_presample_start, quartertodate("1983-Q1"))
s <= Setting(:date_mainsample_start, quartertodate("1983-Q1"))
s <= Setting(:date_mainsample_end, quartertodate("2002-Q4"))
s <= Setting(:date_zlbregime_start, quartertodate("2002-Q4"))

s <= Setting(:n_particles, 5)
s <= Setting(:n_Φ, 10)
s <= Setting(:λ, 2.0)
s <= Setting(:n_smc_blocks, 1)

# Reading in reference values
acc_loglh = h5read("$path/../reference/mutation_RWMH.h5","Acc_loglh")
acc_post = h5read("$path/../reference/mutation_RWMH.h5","Acc_post")

rej_loglh = h5read("$path/../reference/mutation_RWMH.h5","Rej_loglh")
rej_post = h5read("$path/../reference/mutation_RWMH.h5","Rej_post")

prior_draw = h5read("$path/../reference/mutation_RWMH.h5","initial_para")
data = h5read("$path/../reference/smc.h5","data")

out = posterior!(s, prior_draw, data, φ_smc = 1.)
loglh = out[:like]
post = out[:post]
tempering_schedule = vec([.5 .6])
R = eye(n_parameters(s),n_parameters(s))

rvec1 = ones(n_parameters(s),1)
rvec2 = ones(n_parameters(s),1)/100

para_1, loglh_1, post_1, acpt_1 = mutation_RWMH(s,data,prior_draw,loglh,post,2,R,tempering_schedule,rvec=rvec1,rval=0.0)

para_2, loglh_2, post_2, acpt_2 = mutation_RWMH(s,data,prior_draw,loglh,post,2,R,tempering_schedule,rvec=rvec2,rval=1.0)


@test_approx_eq loglh_1 rej_loglh
@test_approx_eq post_1 rej_post
@test_approx_eq acpt_1 0


@test_approx_eq loglh_2 acc_loglh
@test_approx_eq post_2 acc_post
@test_approx_eq acpt_2 1

nothing

