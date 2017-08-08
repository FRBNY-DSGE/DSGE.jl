"""
```
function solve_inefficiency{S<:AbstractFloat}(φ_new::S, φ_old::S, y_t::Array{S, 1}, 
                            p_error::Array{S,2}, EE::Array{S,2}; initialize::Bool=false)
```
Returns the value of the ineffeciency function (average of the normalized weights squared).
"""
function solve_inefficiency{S<:AbstractFloat}(φ_new::S, φ_old::S, y_t::Array{S,1}, 
                            p_error::Array{S,2}, EE::Array{S,2}; initialize::Bool=false)

    n_particles = size(p_error, 2)
    n_obs       = length(y_t)
    w           = zeros(n_particles)
    inv_EE      = inv(EE)
    det_EE      = det(EE)

    # Inefficiency function during initialization
    if initialize
        for i=1:n_particles
            w[i] = ((φ_new/(2*pi))^(n_obs/2) * (det_EE^(-1/2)) * exp(-1/2 * p_error[:,i]' * 
                                                φ_new * inv_EE * p_error[:,i]))[1]
        end
    
    # Inefficiency function during tempering steps
    else
        for i=1:n_particles
            w[i] = (φ_new/φ_old)^(n_obs/2) * exp(-1/2 * p_error[:,i]' * 
                                                 (φ_new-φ_old) * inv_EE * p_error[:,i])[1]
        end
    end
    W = w/mean(w)
    return sum(W.^2)/n_particles
end