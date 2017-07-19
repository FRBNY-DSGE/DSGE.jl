using DSGE
using HDF5, Base.Test
using DataFrames

path = dirname(@__FILE__)

s = AnSchorfheide(testing=true)

s <= Setting(:n_particles, 5)
s <= Setting(:n_Φ, 10)
s <= Setting(:λ, 2.0)

data = h5read("$path/../../test/reference/mutation_RWMH.h5","data")
para_init = h5read("$path/../../test/reference/mutation_RWMH.h5","prior_draw")
post_init = posterior!(s, para_init, data, φ_smc = 1.)
like_init = post_init - prior(s)
tempering_schedule = vec([.5 .6])
R = eye(n_parameters(s),n_parameters(s))
i = 2

#Test 1
#check base case if block size 1 performs the same as standard mutation_RWMH
s <= Setting(:n_smc_blocks, 1)
s <= Setting(:n_MH_steps_smc, 1)

srand(123)
rvec = randn(n_parameters(s),1)
rval = rand()
b_out1, b_out2, b_out3, b_out4 = mutation_block_RWMH(s,data,para_init,like_init,post_init,i,R,tempering_schedule;rvec=rvec,rval=rval)

srand(123)
rvec = randn(n_parameters(s),1)
rval = rand()
out1, out2, out3, out4 = mutation_RWMH(s,data,para_init,like_init,post_init,i,R,tempering_schedule;rvec=rvec,rval=rval)

@test_approx_eq b_out1 out1
@test_approx_eq b_out2 out2
@test_approx_eq b_out3 out3
@test_approx_eq b_out4 out4

#Test 2
#check if each blocking works properly
s <= Setting(:n_smc_blocks, 5)
s <= Setting(:n_MH_steps_smc, 1)
s <= Setting(:step_size_smc, 1)
rstep = ones(n_parameters(s))
out1, _1, _2, _3 = mutation_block_RWMH(s,data,para_init,like_init,post_init,i,R,tempering_schedule;rstep=rstep)

@test_approx_eq out1-para_init [1.,1.,1.,1.,2.,2.,2.,3.,3.,3.,4.,4.,4.,0.,0.,0.] 

nothing
