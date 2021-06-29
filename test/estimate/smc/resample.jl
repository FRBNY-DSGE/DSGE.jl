writing_output = false
if VERSION < v"1.5"
    ver = "111"
elseif VERSION < v"1.6"
    ver = "150"
else
    ver = "160"
end

@everywhere Random.seed!(42)

weights = rand(400)
weights = weights ./ sum(weights)

test_sys_resample    = SMC.resample(weights, method = :systematic)
test_multi_resample  = SMC.resample(weights, method = :multinomial)
test_poly_resample   = SMC.resample(weights, method = :polyalgo)

if writing_output
    jldopen("reference/resample_version=" * ver * ".jld2",
            true, true, true, IOStream) do file
        write(file, "sys", test_sys_resample)
        write(file, "multi", test_multi_resample)
        write(file, "poly", test_poly_resample)
    end
end

saved_sys_resample   = load("reference/resample_version=" * ver * ".jld2", "sys")
saved_multi_resample = load("reference/resample_version=" * ver * ".jld2", "multi")
saved_poly_resample  = load("reference/resample_version=" * ver * ".jld2", "poly")

####################################################################

@testset "Resampling methods" begin
    @test test_sys_resample   == saved_sys_resample
    @test test_multi_resample == saved_multi_resample
    @test test_poly_resample  == saved_poly_resample
end
