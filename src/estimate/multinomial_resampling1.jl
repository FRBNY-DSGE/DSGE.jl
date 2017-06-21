using StatsBase

function multinomial_resampling1(m,weight)
    numParticles=length(weight)
    idx=sample(collect(1:numParticles), weights(weight), numParticles, replace=true)
    return idx
end

println(multinomial_resampling1(1, [.2,.2,.3,.3, .1]))