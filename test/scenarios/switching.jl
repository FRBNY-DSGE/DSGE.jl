using DSGE, Base.Test

# Test switch
original = zeros(1,4)
default  = ones(1,4)

@test DSGE.switch(original, default, zeros(4), zeros(4)) == default
@test DSGE.switch(original, default, ones(4), zeros(4)) == original
@test DSGE.switch(original, default, ones(4), ones(4)) == [0. 1. 1. 1.]

# Test choose_last_period
@test DSGE.choose_last_period(ones(3), 1, 3) == 0
@test DSGE.choose_last_period(ones(3), 2, 3) == 1
@test DSGE.choose_last_period(ones(3), 4, 3) == 3
