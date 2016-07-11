function multinomial_resampling(w)

np = length(w)

w = w'
cw = cumsum(w)

uu = rand(length(w), 1)
indx = zeros(np, 1)
for i = 1:np
    
    u = uu(i)

    j=1
    while j <= np
       if (u < cw(j)), break, end
       
       j = j+1
        
    end

    indx(i) = j

end

m = 0

return indx, m
end
