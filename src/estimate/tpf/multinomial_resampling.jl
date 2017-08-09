using StatsBase
"""
```
multinomial_resampling{S<:Float64}(particle_weights::Array{S,1})
```
Performs multinomial resampling from vector of particle_weights, returning indices of resampled particles.

### Inputs
- `particle_weights`: vector of weights

### Outputs
- `ids`: indices of resampled particles
"""
function multinomial_resampling{S<:Float64}(particle_weights::Array{S,1})
    # Store length of weight vector
    n_particles = size(particle_weights,1)
    # Resample
    ids = sample(collect(1:n_particles), weights(particle_weights), n_particles, replace=true)
    return ids
end

