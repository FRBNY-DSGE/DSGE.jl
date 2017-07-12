"""
```
ineff_func(φ_new::Float64, φ_old::Float64, yt::Array{Float64, 1}, perror::Array, H::Array{Float64}; initialize::Int64=0)
```
Returns the value of the ineffeciency function (average of the normalized weights squared).
"""

function ineff_func(φ_new::Float64, φ_old::Float64, yt::Array{Float64, 1}, perror::Array{Float64,2}, H::Array{Float64,2}; initialize::Int64=0)
    
    n_particles = size(perror, 2)
    w = zeros(n_particles)
    inv_H = inv(H)
    det_H = det(H)

    # Inefficiency function during tempering steps
    if initialize==0
        for i=1:n_particles
            w[i] = (φ_new/φ_old)^(length(yt)/2)*exp(-1/2*perror[:,i]'*(φ_new-φ_old)*inv_H*perror[:,i])[1]
            #w[i] = density(φ_new, φ_old, yt, perror[:,i], H, initialize=initialize)
        end

    # Inefficiency function during initialization
    else
        for i=1:n_particles
            w[i] = ((φ_new/(2*pi))^(length(yt)/2)*(det_H^(-1/2))*exp((-1/2)*perror[:,i]'*φ_new*inv_H*perror[:,i]))[1]
            #w[i] = density(φ_new, 1.0, yt, perror[:,i], H, initialize=initialize)
        end
    end
    W = w/mean(w)
    return sum(W.^2)/n_particles
end