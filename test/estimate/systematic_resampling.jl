using DSGE

# Set parameters for testing
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")))
m = AnSchorfheide(custom_settings = custom_settings, testing = true)
m<=Setting(:use_parallel_workers,true)

# Set number of draws
draws = 10000

# Remove previous testing file
rm("systematic_resampling_test.csv")

# Append weights vectors for testing
weights_vec=Any[]
append!(weights_vec, [collect(0.1:0.1:0.9)])
append!(weights_vec,[[0.4, 0.1, 0.1, 0.2, 0.3]])
append!(weights_vec,[[0.1, 0.8, 0.1]])

for weights in weights_vec
    
    n_parts = length(weights)

    count = zeros(n_parts)
    out = zeros(n_parts)

    for i=1:draws
        out = systematic_resampling(m, weights)
        for idx in out
            count[idx] += 1
        end
    end

    count = count./(draws*n_parts)
    count = round(count,3)
    act_weights = round(weights./sum(weights),3)

    open("systematic_resampling_test.csv","a") do file
        write(file,"Actual Probabilities \n")
        writecsv(file,act_weights')
        write(file,"\nTested Probabilities \n")
        writecsv(file,count')
        write(file,"\n---------------------\n")
    end
end
nothing