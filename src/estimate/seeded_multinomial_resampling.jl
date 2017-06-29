function seeded_multinomial_resampling(weight::Array{Float64,1})
    num_particles=size(weight,1)
    weight = weight./sum(weight)
    weight *= num_particles
    weight = round(weight)    
    sorted_idx = sortperm(weight,rev=true)
    out_weight::Array{Int64, 1}=zeros(num_particles)
    i=k=is_first=1

    while i<=num_particles
        if (weight[sorted_idx[k]]>0.0) | (is_first==1)
            weight[sorted_idx[k]]-=1
            out_weight[i]=Integer(sorted_idx[k])
            i+=1
            is_first=0
        else
            k+=1
            is_first=1
        end
    end
    return out_weight
end

