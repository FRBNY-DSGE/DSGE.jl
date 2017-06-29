using StatsBase

function multinomial_resampling(weight::Array{Float64,1})
    num_particles=size(weight,1)
    ids=sample(collect(1:num_particles), weights(weight), num_particles, replace=true)
    return ids
end

