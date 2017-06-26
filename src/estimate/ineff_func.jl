function ineff(φn::Float64, φn_1::Float64, yt::Array{Float64, 1}, perror::Array{Float64}, H::Array{Float64}; initialize::Int64=0)
    N = size(perror, 2)
    w=Array{Float64,1}
    # Inefficiency function during tempering steps
    if initialize==0
        for i=1:N
            w[i] = density(φn, φn_1, yt, perror[:,i], H, initialize=1)
        end
    #Inefficiency function during initialization (note φn_1 is never called)
    else
        for i=1:N
            w[i] = density(φn, 1, yt, perror[:,i], H, initialize=0)
        end
    end
    W = w/mean(w)
    return sum(W.^2)/N
end