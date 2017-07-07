"""
```
systematic_resampling(m, weight)
```

Reindexing and reweighting samples from a degenerate distribution

### Arguments:
- `weight`: wtsim[:,i]
        the weights of a degenerate distribution.

### Output:

- `vec(indx)`: id
        the newly assigned indices of parameter draws.

"""
function systematic_resampling(m::AbstractModel, weight::AbstractArray)
    
    n_part = length(weight)
    weight = weight./sum(weight)

    # Stores cumulative weights until given index
    cumulative_weight = cumsum(weight)
    weight = weight'
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
            if (u < cumulative_weight[j])
                break
            end
            j = j+1
        end
        indx[i] = j
    end

    # Map function if parallel
    parallel = get_setting(m,:use_parallel_workers)
    if parallel
        parindx = pmap(j -> subsys(j), 1:n_part)
    else
        parindx = [subsys(j) for j = 1:n_part]'
    end

    # Transpose and round output indices
    indx = parindx'
    indx = round(Int, indx)
    
    # Write output to file if in testing mode
    if m.testing
        open("resamples.csv","a") do x
            writecsv(x,indx')
        end
    end
    return vec(indx)
end
