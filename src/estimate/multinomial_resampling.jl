using StatsBase

function multinomial_resampling(weight::Array{Float64,1})
    numParticles=size(weight,1)
    idx=sample(collect(1:numParticles), weights(weight), numParticles, replace=true)
    return idx
end

