using DSGE, Base.Test

# Test switch
original = zeros(1,4)
default  = ones(1,4)

@test DSGE.switch(original, default, zeros(4), zeros(4)) == default
@test DSGE.switch(original, default, ones(4), zeros(4)) == original
@test DSGE.switch(original, default, ones(4), ones(4)) == [0. 1. 1. 1.]

# Test choose_last_period
@test DSGE.choose_switch_period(ones(3), 1) == 1
@test DSGE.choose_switch_period(ones(3), 2) == 2
@test DSGE.choose_switch_period(ones(3), 5) == 4
@test DSGE.choose_switch_period(zeros(3), 1) == 4
