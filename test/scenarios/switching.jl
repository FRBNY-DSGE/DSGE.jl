using DSGE, Test

# Test switch
original = zeros(1,4)
default  = ones(1,4)

@testset "Test scenario switch" begin
    @test DSGE.switch(original, default, zeros(4), zeros(4))[1] == default
    @test DSGE.switch(original, default, ones(4), zeros(4))[1] == original
    @test DSGE.switch(original, default, ones(4), ones(4))[1] == [0. 1. 1. 1.]
end

# Test choose_last_period
@testset "Test scenario switch in the last period" begin
    @test DSGE.choose_switch_period(ones(3), 1) == 1
    @test DSGE.choose_switch_period(ones(3), 2) == 2
    @test DSGE.choose_switch_period(ones(3), 5) == 4
    @test DSGE.choose_switch_period(zeros(3), 1) == 4
end
