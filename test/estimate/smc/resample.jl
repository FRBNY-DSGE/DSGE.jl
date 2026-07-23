using BenchmarkTools

writing_output = false
# RNG-dependent references: Julia 1.7+ switched the default RNG to a per-Task Xoshiro256++,
# so the random weights and seeded resampling draws differ from the "150"/"160" data —
# regenerate with writing_output on under the target Julia.
if VERSION < v"1.5"
    ver = "111"
elseif VERSION < v"1.6"
    ver = "150"
elseif VERSION < v"1.7"
    ver = "160"
else
    ver = "1126"
end

path = dirname(@__FILE__)

# Plain seed (not @everywhere): in a single process @everywhere doesn't pin the task-local
# RNG the seeded draws use, leaving the RNG references unreproducible.
Random.seed!(42)

weights = rand(400)
weights = weights ./ sum(weights)

test_sys_resample    = SMC.resample(weights, method = :systematic)
test_multi_resample  = SMC.resample(weights, method = :multinomial)
test_poly_resample   = SMC.resample(weights, method = :polyalgo)

if writing_output
    JLD2.jldopen("$path/../../reference/resample_version=" * ver * ".jld2",
            true, true, true, IOStream) do file
        write(file, "sys", test_sys_resample)
        write(file, "multi", test_multi_resample)
        write(file, "poly", test_poly_resample)
    end
end

saved_sys_resample   = load("$path/../../reference/resample_version=" * ver * ".jld2", "sys")
saved_multi_resample = load("$path/../../reference/resample_version=" * ver * ".jld2", "multi")
saved_poly_resample  = load("$path/../../reference/resample_version=" * ver * ".jld2", "poly")

####################################################################

@testset "Resampling methods" begin
    @test test_sys_resample   == saved_sys_resample
    @test test_multi_resample == saved_multi_resample
    @test test_poly_resample  == saved_poly_resample
end

################
# Benchmarking #
################
# Flip to true to run; off by default. Pure numerics, no FRED API.
run_benchmarks = false

if run_benchmarks
    b_sys   = @benchmark SMC.resample($weights, method = :systematic)
    b_multi = @benchmark SMC.resample($weights, method = :multinomial)
    b_poly  = @benchmark SMC.resample($weights, method = :polyalgo)

    println("\n===== estimate/smc/resample benchmark results =====")
    for (name, b) in [("resample systematic  (n=400)", b_sys),
                      ("resample multinomial (n=400)", b_multi),
                      ("resample polyalgo    (n=400)", b_poly)]
        println(name, "  time:   ", BenchmarkTools.prettytime(median(b).time),
                "   memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end

nothing
