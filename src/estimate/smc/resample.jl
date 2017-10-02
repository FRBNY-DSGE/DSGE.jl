"""
```
resample(weights::AbstractArray; method::Symbol=:systematic,
         parallel::Bool=false, testing::Bool=false)
```

Reindexing and reweighting samples from a degenerate distribution

### Arguments:
- `weights`: get_weights(cloud)
        the weights of a degenerate distribution.
- `method`: :systematic or :multinomial
        the method for resampling
- `parallel`: to indicate whether to resample using multiple workers (if available)
- `testing`: to indicate whether to give test output

### Output:

- `vec(indx)`: id
        the newly assigned indices of parameter draws.
"""
function resample(weights::AbstractArray; method::Symbol=:systematic,
                  parallel::Bool=false, testing::Bool=false)
    if method == :systematic
        n_parts = length(weights)
        # Stores cumulative weights until given index
        cumulative_weights = cumsum(weights./sum(weights))
        uu = Vector{Float64}(n_parts)

        # Random part of algorithm - choose offset of first index by some u~U[0,1)
        rand_offset = rand()

        # Set "spokes" at the position of the random offset
        for j = 1:n_parts
            uu[j] = (j - 1) + rand_offset
        end

        # Function solves where an individual "spoke" lands
        function subsys(i::Int)
            findfirst(j -> cumulative_weights[j] > uu[i]/n_parts, 1:length(cumulative_weights))
        end

        # Map function if parallel
        if parallel
            indx =
            @sync @parallel (vcat) for j in 1:n_parts
                subsys(j)
            end
        else
            indx = [subsys(j) for j = 1:n_parts]
        end

        return indx

    elseif method == :multinomial
        n_part = length(weights)
        weights = Weights(weights./sum(weights))
        # to use the proper method of sample, weight needs to be an AbstractWeight object
        return sample(1:n_part,weights,n_part,replace=true)
    else
        throw("Invalid resampler in SMC. Set model setting :resampler_smc to either :systematic or :multinomial")
    end
end
