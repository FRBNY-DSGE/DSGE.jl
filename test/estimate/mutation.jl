using DSGE
using JLD, HDF5, Base.Test
using DataFrames

path = dirname(@__FILE__)
# path = "/home/rcemxc29/.julia/v0.5/DSGE/test/estimate"

custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start => Setting(:date_forecast_start, quartertodate("2015-Q4")))
s = AnSchorfheide(custom_settings = custom_settings, testing = true)

s <= Setting(:n_particles, 5)
s <= Setting(:n_Φ, 10)
s <= Setting(:λ, 2.0)
cloud = ParticleCloud(s, get_setting(s, :n_particles))
tempering_schedule = [.5,.6]

file = jldopen("$path/../reference/mutation.jld","r")
# file = jldopen("mutation.jld","r")
prev_part = read(file, "prev_part")
good_part = read(file, "good_part")
bad_part  = read(file, "bad_part")
rvec1     = read(file, "rvec1")
rvec2     = read(file, "rvec2")
data      = read(file, "data")
close(file)

# Test 1
current_φ = tempering_schedule[2]
previous_φ = tempering_schedule[1]
R = eye(n_parameters(s))

# bad particle
test_prev_part = deepcopy(prev_part)
update!(s, test_prev_part.value)
prev_part.logpost = previous_φ*prev_part.loglh + prior(s)
test_prev_part.logpost = previous_φ*test_prev_part.loglh + prior(s)
new_part_0 = mutation(s,data,test_prev_part,R,current_φ,previous_φ,rvec=rvec2,rval=1.0)

@test_approx_eq test_prev_part.value prev_part.value #testing if prev_part changes in place due to update_mutation!
@test_approx_eq new_part_0.loglh prev_part.loglh
@test_approx_eq new_part_0.logpost (prev_part.logpost + (current_φ - previous_φ)*prev_part.loglh)
@test_approx_eq new_part_0.accept false

# good particle
test_prev_part = deepcopy(prev_part)
update!(s, good_part.value)
prev_part.logpost = previous_φ*prev_part.loglh + prior(s)
good_part.logpost = current_φ*good_part.loglh + prior(s)
new_part_1 = mutation(s,data,test_prev_part,R,current_φ,previous_φ,rvec=rvec1,rval=0.0)

@test_approx_eq new_part_1.loglh good_part.loglh
@test_approx_eq new_part_1.logpost good_part.logpost
@test_approx_eq new_part_1.accept true

# Test 2
# check if each blocking works properly
s <= Setting(:n_smc_blocks, 5)
s <= Setting(:n_MH_steps_smc, 1)
s <= Setting(:step_size_smc, 1.)
fixed_inds = [θ.fixed for θ in s.parameters]
rstep = ones(n_parameters(s))
current_φ = tempering_schedule[2]
previous_φ = tempering_schedule[1]
R = eye(n_parameters(s))
test_prev_part = deepcopy(prev_part)
new_part = mutation(s, data, test_prev_part, R, current_φ, previous_φ; rstep = rstep)

@test_approx_eq new_part.value-prev_part.value [1.,1.,1.,1.,2.,2.,2.,3.,3.,3.,4.,4.,4.,5.,5.,5.]

nothing
