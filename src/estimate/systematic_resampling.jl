function systematic_resampling( w )

np = length(w)
w = w'
cw = cumsum(w)
uu = zeros(length(w),1)
csi=rand(1)

for j=1:length(w)
    uu(j) = (j-1)+csi
end
    
indx = zeros(np, 1)

@parallel for i = 1:np
    u = uu(i)/length(w)
    j=1
    while j <= np
        if (u < cw(j)) 
            break
        end
        j = j+1
    end
    indx(i) = j
end

m = 0

return indx, m
end
