"""
```
systematic_resampling(wght)
```

Reindexing and reweighting samples from a degenerate distribution

### Arguments:
- `wght`: wtsim[:,i]
        the weights of a degenerate distribution.

### Output:

- `vec(indx)`: id
        the newly assigned indices of parameter draws.

"""
function systematic_resampling(wght)


npart = length(wght)
wght = wght'
cwght = cumsum(wght')
uu = zeros(npart,1)
csi=rand()

for j=1:npart
    uu[j] = (j-1)+csi
end
    
indx = zeros(npart, 1)

#@parallel 
#for i = 1:npart
#    u = uu[i]/npart
#    j=1
#    while j <= npart
#        if (u < cwght[j]) 
#            break
#        end
#        j = j+1
#    end
#    indx[i] = j
#end

function subsys(i)
    u = uu[i]/npart
    j=1
    while j <= npart
        if (u < cwght[j]) 
            break
        end
        j = j+1
    end
    indx[i] = j
end


parindx = @sync @parallel (hcat) for j = 1:npart 
    subsys(j)
end
indx = parindx'

indx = round(Int, indx)

open("resamples.csv","a") do x
    writecsv(x,indx')
end 


return vec(indx)
end
