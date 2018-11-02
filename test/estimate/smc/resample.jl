# To be removed after running this test individually in the REPL successfully
using DSGE
using HDF5, JLD2
using Test

Random.seed!(42)
weights = rand(400)
weights = weights/sum(weights)
test_sys_resample = resample(weights, method = :systematic)
test_multi_resample = resample(weights, method = :multinomial)
test_poly_resample = resample(weights, method = :polyalgo)

saved_sys_resample = load("reference/resample.jld2", "sys")
saved_multi_resample = load("reference/resample.jld2", "multi")
saved_poly_resample = load("reference/resample.jld2", "poly")

####################################################################
@testset "Resampling methods" begin
    @test test_sys_resample == saved_sys_resample
    @test test_multi_resample == saved_multi_resample
    @test test_poly_resample == saved_poly_resample
end
