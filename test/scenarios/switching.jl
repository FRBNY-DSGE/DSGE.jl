using DSGE, JLD

# Test switch()
original = zeros(1,4)
default  = ones(1,4)

@tests switch(original, default, zeros(4), zeros(4)) == default
@tests switch(original, default, ones(4), zeros(4)) == original
@tests switch(original, default, ones(4), zeros(4)) == [0. 1. 1. 1. 1.]

# Test choose_last_period
@tests choose_last_period(ones(3), 1, 3) == 0
@tests choose_last_period(ones(3), 2, 3) == 1
@tests choose_last_period(ones(3), 4, 3) == 3
