"""
```
solve_inefficiency{S<:AbstractFloat}(φ_new::S, φ_old::S, y_t::Array{S, 1}, p_error::Array{S,2}, 
                                     HH::Array{S,2}; initialize::Bool=false)
```
Returns the value of the ineffeciency function InEff(φₙ), where:

        InEff(φₙ) = (1/M) ∑ᴹ (W̃ₜʲ(φₙ))² 
                    
Where ∑ is over j=1...M particles, and for a particle j: 

        W̃ₜʲ(φₙ) = w̃ₜʲ(φₙ) / (1/M) ∑ᴹ w̃ₜʲ(φₙ)

Where ∑ is over j=1...M particles, and incremental weight is:
            
        w̃ₜʲ(φₙ) = pₙ(yₜ|sₜʲ'ⁿ⁻¹) / pₙ₋₁(yₜ|sₜ^{j,n-1})
                = (φₙ/φₙ₋₁)^(d/2) exp{-1/2 [yₜ-Ψ(sₜʲ'ⁿ⁻¹)]' (φₙ-φₙ₋₁) ∑ᵤ⁻¹ [yₜ-Ψ(sₜʲ'ⁿ⁻¹)]}

### Inputs

- `φ_new`: φₙ
- `φ_old`: φₙ₋₁
- `y_t`: vector of observables for time t
- `p_error`: (`n_states` x `n_particles`) matrix of particles' errors yₜ - Ψ(sₜʲ'ⁿ⁻¹) in columns
- `HH`: measurement error covariance matrix, ∑ᵤ

### Keyword Arguments

- `initialize::Bool`: flag to indicate whether this is being used in initialization stage, 
    in which case one instead solves the formula for w̃ₜʲ(φₙ) as:

    w̃ₜʲ(φ₁) = (φ₁/2π)^(d/2)|∑ᵤ|^(1/2) exp{-1/2 [yₜ-Ψ(sₜʲ'ⁿ⁻¹)]' φ₁ ∑ᵤ⁻¹ [yₜ-Ψ(sₜʲ'ⁿ⁻¹)]}

"""
function solve_inefficiency_refactor{S<:AbstractFloat}(φ_new::S, φ_old::S, y_t::Array{S,1}, 
                            p_error::Array{S,2}, HH::Array{S,2}; initialize::Bool=false)

    n_particles = size(p_error, 2)
    n_obs       = length(y_t)
    w           = zeros(n_particles)
    inv_HH      = inv(HH)
    det_HH      = det(HH)

    # Inefficiency function during initialization
    if initialize
    top=zeros(n_particles)
    bottom=zeros(n_particles)
    for i=1:n_particles
#changed from φ_new/2*pi to φ_new/φ_old and got rid of det_HH⁽-1/2) and made φ_new-φ_old from just φ_new
            w[i] = ((φ_new)/(φ_old))^(n_obs/2) * (exp(-1/2 * p_error[:,i]' * 
                                                (φ_new-φ_old) * inv_HH * p_error[:,i]))[1]
            top[i] = exp(-2.0*(φ_new-φ_old)*p_error[:,i]'*inv_HH*p_error[:,i])[1]
            bottom[i] = exp(-(φ_new-φ_old)*p_error[:,i]'*inv_HH*p_error[:,i])[1]
        end
    
    # Inefficiency function during tempering steps
    else
        for i=1:n_particles
            w[i] = (φ_new/φ_old)^(n_obs/2) * exp(-1/2 * p_error[:,i]' * 
                                                 (φ_new-φ_old) * inv_HH * p_error[:,i])[1]
        end
    end

   # return mean(top)/(mean(bottom)^2)
    W = w/mean(w)
    return sum(W.^2)/n_particles
end