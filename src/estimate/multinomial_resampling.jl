using StatsBase

function multinomial_resampling1(m::AbstractModel,weight::Array{Float64,1})
    numParticles=length(weight)
    idx=sample(collect(1:numParticles), weights(weight), numParticles, replace=true)
    return idx
end

