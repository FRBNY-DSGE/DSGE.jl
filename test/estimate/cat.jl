using BenchmarkTools

path = dirname(@__FILE__)

# Set up
m = AnSchorfheide(testing = true)

kal1, kal2 = JLD2.jldopen("$path/../reference/kalman_cat_args.jld2", "r") do file
    read(file, "kal1"), read(file, "kal2")
end

# Concatenate Kalmans
kal12 = cat(m, kal1, kal2)

# Test equality
exp_kal12 = JLD2.jldopen("$path/../reference/kalman_cat_out.jld2", "r") do file
    read(file, "kal12")
end

@testset "Testing Kalman output concatenation" begin
    for arg in fieldnames(typeof(kal1))
        @test exp_kal12[arg] ≈ kal12[arg]
    end
end

################
# Benchmarking #
################
# Flip to true to run; off by default. Inputs are local JLD2, no FRED API.
run_benchmarks = false

if run_benchmarks
    b_cat = @benchmark cat($m, $kal1, $kal2)

    println("\n===== estimate/cat benchmark results =====")
    println("cat (Kalman)  time:   ", BenchmarkTools.prettytime(median(b_cat).time),
            "   memory: ", BenchmarkTools.prettymemory(median(b_cat).memory))
end

nothing
