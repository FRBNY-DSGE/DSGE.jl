using StatsBase
function multinomial_sampling2(m::AbstractModel,w)
	 n = length(w)
	 idx = zeros(n, 1)
	 weights = StatsBase.ProbabilityWeights(w)	  
	 for i = 1:n
	     idx[i] = StatsBase.sample(collect(1:n),weights)
	 end
	 return idx
end