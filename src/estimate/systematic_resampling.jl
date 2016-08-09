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

m = 0

indx = round(Int, indx)
return vec(indx), m
end
