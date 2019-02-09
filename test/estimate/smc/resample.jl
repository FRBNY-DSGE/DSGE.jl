# To be removed after running this test individually in the REPL successfully
using DSGE
using HDF5, JLD
import Base.Test: @test, @testset

@everywhere srand(42)
weights = rand(400)
weights = weights/sum(weights)
test_sys_resample = resample(weights, method = :systematic)
test_multi_resample = resample(weights, method = :multinomial)
test_poly_resample = resample(weights, method = :polyalgo)

#=jldopen("resample.jld2", "w") do file
    file["sys"] = test_sys_resample
    file["multi"] = test_multi_resample
    file["poly"] = test_poly_resample
end=#

saved_sys_resample = load("reference/resample.jld2", "sys")
saved_multi_resample = load("reference/resample.jld2", "multi")
saved_poly_resample = load("reference/resample.jld2", "poly")

####################################################################
@testset "Resampling methods" begin
    @test test_sys_resample == saved_sys_resample
    @test test_multi_resample == saved_multi_resample
    @test test_poly_resample == saved_poly_resample
end
