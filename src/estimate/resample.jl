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
        n_part = length(weights)
        weights = weights./sum(weights)
        # Stores cumulative weights until given index
        cumulative_weights = cumsum(weights)
        weights = weights'
        uu = zeros(n_part,1)

        # Random part of algorithm - choose offset of first index by some u~U[0,1)
        rand_offset=rand()

        # Set "spokes" at the position of the random offset
        for j=1:n_part
            uu[j] = (j-1)+rand_offset
        end

        # Initialize output vector
        indx = zeros(n_part, 1)

        # Function solves where an individual "spoke" lands
        function subsys(i)
            u = uu[i]/n_part
            j=1
            while j <= n_part
                if (u < cumulative_weights[j])
                    break
                end
                j = j+1
            end
            indx[i] = j
        end

        # Map function if parallel
        if parallel
            parindx =
            @parallel (vcat) for j in 1:n_part
                subsys(j)
            end
        else
            parindx = [subsys(j) for j = 1:n_part]'
        end

        # Transpose and round output indices
        indx = parindx'
        indx = round(Int, indx)

        return vec(indx)
    elseif method == :multinomial
        n_part = length(weights)
        weights = Weights(weights./sum(weights))
        # to use the proper method of sample, weight needs to be an AbstractWeight object
        return sample(1:n_part,weights,n_part,replace=true)
    else
        throw("Invalid resampler in SMC. Set model setting :resampler_smc to either :systematic or :multinomial")
    end
end
