"""
```
ineff_func(φ_new::Float64, φ_old::Float64, y_t::Array{Float64, 1}, perror::Array, EE::Array{Float64}; initialize::Int64=0)
```
Returns the value of the ineffeciency function (average of the normalized weights squared).
"""

function ineff_func(φ_new::Float64, φ_old::Float64, y_t::Array{Float64, 1}, perror::Array{Float64,2}, EE::Array{Float64,2}; initialize::Int64=0)
    
    n_particles = size(perror, 2)
    n_obs       = length(y_t)
    w           = zeros(n_particles)
    inv_EE      = inv(EE)
    det_EE      = det(EE)

    # Inefficiency function during tempering steps
    if initialize==0
        for i=1:n_particles
            w[i] = (φ_new/φ_old)^(n_obs/2)*exp(-1/2*perror[:,i]'*(φ_new-φ_old)*inv_EE*perror[:,i])[1]
        end

    # Inefficiency function during initialization
    else
        for i=1:n_particles
            w[i] = ((φ_new/(2*pi))^(n_obs/2)*(det_EE^(-1/2))*exp((-1/2)*perror[:,i]'*φ_new*inv_EE*perror[:,i]))[1]
        end
    end
    W = w/mean(w)
    return sum(W.^2)/n_particles
end