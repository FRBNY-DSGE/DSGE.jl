"""
```
resample(weights::AbstractArray; method::Symbol = :systematic)
```

Reindexing and reweighting samples from a degenerate distribution

### Arguments:
- `weight`: wtsim[:,i]
        the weights of a degenerate distribution.
- `method`: :systematic, :multinomial, or :polyalgo
        the method for resampling

### Output:

- `indx`: the newly assigned indices of parameter draws.
"""
function resample(weights::Vector{Float64}; method::Symbol = :systematic)
    if method == :multinomial
        n_parts = length(weights)
        indx = Vector{Int64}(n_parts)

        # Stores cumulative weights until given index
        cumulative_weights = cumsum(weights/sum(weights))
        offset = rand(n_parts)

        # This could be parallelized
        for i in 1:n_parts
            indx[i] = findfirst(x -> offset[i] < x, cumulative_weights)
        end

        return indx
    elseif method == :systematic
        n_parts = length(weights)

        # Stores cumulative weights until given index
        cumulative_weights = cumsum(weights/sum(weights))
        offset = rand()

        # Function solves where an individual "spoke" lands
        function subsys(i::Int, offset::Float64, n_parts::Int64, start_ind::Int64,
                        cumulative_weights::Vector{Float64})
            threshold = (i - 1 + offset)/n_parts
            range = start_ind:n_parts
            for j in range
                if cumulative_weights[j] > threshold
                    return j
                end
            end
            return 0
        end

        indx = Vector{Int64}(n_parts)
        for i in 1:n_parts
            if i == 1
                indx[i] = subsys(i, offset, n_parts, 1, cumulative_weights)
            else
                indx[i] = subsys(i, offset, n_parts, indx[i-1], cumulative_weights)
            end
        end

        return indx
    elseif method == :polyalgo
        n_parts = length(weights)
        weights = Weights(weights./sum(weights))
        return sample(1:n_parts, weights, n_parts, replace = true)
    else
        throw("Invalid resampler in SMC. Set model setting :resampler_smc to either :systematic, :multinomial, or :polyalgo")
    end
end
