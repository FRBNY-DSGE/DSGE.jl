"""
```
ineff_func(φ_new::Float64, φ_old::Float64, yt::Array{Float64, 1}, perror::Array, H::Array{Float64}; initialize::Int64=0)
```
Returns the value of the ineffeciency function (average of the normalized weights squared).
"""

function ineff_func(φ_new::Float64, φ_old::Float64, yt::Array{Float64, 1}, perror::Array, H::Array{Float64}; initialize::Int64=0)
    N = size(perror, 2)
    w=zeros(N)
    # Inefficiency function during tempering steps
    if initialize==0
        for i=1:N
            w[i] = density(φ_new, φ_old, yt, perror[:,i], H, initialize=initialize)
        end
    #Inefficiency function during initialization (note φn_1 is never called)
    else
        for i=1:N
           w[i] = density(φ_new, 1.0, yt, perror[:,i], H, initialize=initialize)
        end
    end
    W = w/mean(w)
    return sum(W.^2)/N
end