
using DSGE
using HDF5, Base.Test
using DataFrames

path = dirname(@__FILE__)

# Set up models for testing
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")))
s = AnSchorfheide(custom_settings = custom_settings, testing = true)

s <= Setting(:n_particles, 5)
s <= Setting(:n_Φ, 10)
s <= Setting(:λ, 2.0)
s <= Setting(:n_smc_blocks, 1)

# Reading in reference values
acc_loglh = h5read("$path/../reference/mutation_RWMH.h5","acc_loglh")
acc_post = h5read("$path/../reference/mutation_RWMH.h5","acc_post")

rej_loglh = h5read("$path/../reference/mutation_RWMH.h5","rej_loglh")
rej_post = h5read("$path/../reference/mutation_RWMH.h5","rej_post")

prior_draw = h5read("$path/../reference/mutation_RWMH.h5","prior_draw")
data = h5read("$path/../reference/mutation_RWMH.h5","data")

post = posterior!(s, prior_draw, data, φ_smc = 1.)
loglh = post - prior(s)
tempering_schedule = vec([.5 .6])
R = eye(n_parameters(s),n_parameters(s))

rvec1 = ones(n_parameters(s),1)
rvec2 = ones(n_parameters(s),1)/100

para_1, loglh_1, post_1, acpt_1 = mutation_RWMH(s,data,prior_draw,loglh,post,2,R,tempering_schedule,rvec=rvec1,rval=0.0)

para_2, loglh_2, post_2, acpt_2 = mutation_RWMH(s,data,prior_draw,loglh,post,2,R,tempering_schedule,rvec=rvec2,rval=1.0)

para_3, loglh_3, post_3, acpt_3 = mutation_block_RWMH(s,data,prior_draw,loglh,post,2,R,tempering_schedule,rvec=rvec2,rval=1.0)


@test_approx_eq loglh_1 rej_loglh
@test_approx_eq post_1 rej_post
@test_approx_eq acpt_1 0


@test_approx_eq loglh_2 acc_loglh
@test_approx_eq post_2 acc_post
@test_approx_eq acpt_2 1

nothing

