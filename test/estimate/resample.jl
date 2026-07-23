using DSGE, DelimitedFiles, BenchmarkTools

path = dirname(@__FILE__)
# Set parameters for testing
custom_settings = [DSGE.Setting(:date_forecast_start, quartertodate("2015-Q4"))]
m = AnSchorfheide(custom_settings = custom_settings, testing = true)
m <= DSGE.Setting(:use_parallel_workers, true)
parallel = DSGE.get_setting(m, :use_parallel_workers)

# Set number of draws
draws = 10000

# Remove previous testing file
rm("$path/../reference/systematic_resampling_test.csv")

# Append weights vectors for testing
weights_vec=Any[]
# Increasing weight
append!(weights_vec, [collect(0.1:0.1:0.9)])
# Even distribution between indices
append!(weights_vec,[[0.4, 0.1, 0.1, 0.2, 0.3]])
# Index in the middle should appear far more frequently than others
append!(weights_vec,[[0.1, 0.8, 0.1]])
# Weights whose sum far exceeds 0
append!(weights_vec,[[9.0, 20.0, 2.0, 7.0, 14.0, 23.0, 4.0]])
# Index with 0 probability
append!(weights_vec,[[0.0, 0.5, 0.5]])

open("$path/../reference/systematic_resampling_test.csv","w") do file
   write(file,"---------------------")
end

for weights in weights_vec

    n_parts = length(weights)

    count = zeros(n_parts)
    out = zeros(n_parts)

    for i=1:draws
        out = resample(weights; method = :systematic, parallel = parallel)
        for idx in out
            count[idx] += 1
        end
    end

    count = count./(draws*n_parts)
    count = round.(count; digits=3)
    act_weights = round.(weights./sum(weights); digits=3)

    open("$path/../reference/systematic_resampling_test.csv","a") do file
        write(file,"\nActual Probabilities: ")
        writedlm(file,act_weights',',')
        write(file,"Tested Probabilities: ")
        writedlm(file,count',',')
        write(file,"---------------------")
    end

end

################
# Benchmarking #
################
# Flip to true to run; off by default. Pure numerics, no FRED API.
run_benchmarks = false

if run_benchmarks
    # Particle-filter-scale weight vector (n=1000, matching n_particles).
    w_large = collect(1.0:1000.0)

    b_large = @benchmark resample($w_large; method = :systematic, parallel = false)

    println("\n===== estimate/resample benchmark results =====")
    println("resample systematic (n=1000)  time:   ", BenchmarkTools.prettytime(median(b_large).time),
            "   memory: ", BenchmarkTools.prettymemory(median(b_large).memory))
end

nothing
