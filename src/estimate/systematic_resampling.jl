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
function systematic_resampling(m, weight)


    npart = length(weight)
    weight = weight'
    cweight = cumsum(weight')
    uu = zeros(npart,1)
    csi=rand()

    for j=1:npart
        uu[j] = (j-1)+csi
    end

    indx = zeros(npart, 1)

    function subsys(i)
        u = uu[i]/npart
        j=1
        while j <= npart
            if (u < cweight[j])
                break
            end
            j = j+1
        end
        indx[i] = j
    end

    parallel = get_setting(m,:use_parallel_workers)
    if parallel
        parindx = @sync @parallel (hcat) for j = 1:npart
            subsys(j)
        end
    else
        parindx = [subsys(j) for j = 1:npart]'
    end
    indx = parindx'

    indx = round(Int, indx)

    if m.testing
        open("resamples.csv","a") do x
            writecsv(x,indx')
        end
    end

    return vec(indx)
end
