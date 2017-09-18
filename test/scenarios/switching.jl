using DSGE, JLD

# Test switch()
original = zeros(1,4)
default  = ones(1,4)

@assert switch(original,default,[0. for x in 1:4],[0. for x in 1:4]) == default
@assert switch(original,default,[1. for x in 1:4],[0. for x in 1:4]) == original
@assert switch(original,default,[1. for x in 1:4],[0. for x in 1:4]) == [0. 1. 1. 1. 1.]

# Test choose_last_period
@assert choose_last_period([1. 1. 1.], 1, 3) == 0
@assert choose_last_period([1. 1. 1.], 2, 3) == 1
@assert choose_last_period([1. 1. 1.], 4, 3) == 3
