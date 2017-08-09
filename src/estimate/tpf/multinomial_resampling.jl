using StatsBase
"""
```
function multinomial_resampling{S<:Float64}(weight::Array{S,1})
```
### Inputs
- `weight`: vector of weigts

### Outputs
- `ids`: indices of resampled particles
"""
function multinomial_resampling{S<:Float64}(weight::Array{S,1})
    # Store length of weight vector
    n_particles = size(weight,1)
    # Resample
    ids = sample(collect(1:n_particles), weights(weight), n_particles, replace=true)
    return ids
end

