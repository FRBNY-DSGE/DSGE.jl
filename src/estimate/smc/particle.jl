"""
```
Particle
```

The `Particle` type contains the values and weight of a given vector of parameters.

### Fields
- `weights::Array{Float64,1}`: The weights of a given particle (from the perspective of importance sampling). The incremental weight that updates the weights at each step φ_n is given by p(Y|θ^i_{n-1})^(φ_n-φ_{n-1}).
- `keys::Array{Symbol,1}`: The key of the parameter corresponding to a given value.
- `value::Array{Float64,1}`: The matrix of particles (n_params by n_parts)
"""
type Particle
    weight::Float64
    keys::Array{Symbol,1}
    value::Array{Float64,1}
    loglh::Float64
    logpost::Float64
    accept::Bool
end

"""
```
ParticleCloud
```

The `ParticleCloud` type contains all of the relevant information for a given cloud of particles
in the SMC algorithm. Information for a single iteration is stored at any given time (and thus
the final output will be the final cloud of particles, of which only the particle values will be saved).

### Fields
- `particles::Array{Particle,1}`: The array of particles (which contain weight, keys, and value information)
- `n_Φ::Int`: The stage of the tempering schedule
- `rsmp::Int`: A binary indicator of whether or not the particles were resampled in a given stage
- `ESS::Float64`: The effective sample size (resample if ESS falls under n_parts/2)
- `accept::Float64`: The average acceptance rate of mutation steps
"""
type ParticleCloud
    particles::Array{Particle,1}
    tempering_schedule::Array{Float64,1}
    stage_index::Int
    n_Φ::Int
    resamples::Int
    c::Float64
    ESS::Float64
    accept::Float64
end

# Easier constructor for ParticleCloud, which initializes the weights to be equal, and everything else
# in the Particle object etc. to be empty
function ParticleCloud(m::AbstractModel, n_parts::Int)
    return ParticleCloud([Particle(1/n_parts,[m.parameters[i].key for i in 1:length(m.parameters)],
                         zeros(length(m.parameters)),0.,0.,false) for n in 1:n_parts],
                         zeros(get_setting(m,:n_Φ)),1,0,0,0.,0.,0.)
end

function get_weights(c::ParticleCloud)
    return map(p -> p.weight,c.particles)
end

function get_vals(c::ParticleCloud)
    temp = map(p -> p.value, c.particles)
    mat = temp[1]
    for i in 2:length(temp)
        mat = hcat(mat,temp[i])
    end
    return mat
end

function get_loglh(c::ParticleCloud)
    temp = map(p -> p.loglh, c.particles)
    arr = temp[1]
    for i in 2:length(temp)
        arr = vcat(arr,temp[i])
    end
    return arr
end

function get_logpost(c::ParticleCloud)
    temp = map(p -> p.logpost, c.particles)
    arr = temp[1]
    for i in 2:length(temp)
        arr = vcat(arr,temp[i])
    end
    return arr
end

function update_draws!(c::ParticleCloud, draws::Array{Float64,2})
    # transform matrix draws into array of rows vectors
    # for the purposes of aligning zip on the right dimension
    @assert size(draws) == (length(c.particles[1]),length(c))
    draws = [draws[:,i] for i in 1:size(draws)[2]]
    for (p,draw) in zip(c.particles,draws)
        update_val!(p,draw)
    end
end

function update_draws!(c::ParticleCloud, draws::Array{Particle,1})
    for (i,d) in enumerate(draws)
        c.particles[i] = d
    end
end

function update_weights!(c::ParticleCloud, incweight::Array{Float64,1})
    new_weights = get_weights(c).*incweight
    for (p,w) in zip(c.particles, new_weights)
        update_weight!(p,w)
    end
end

function update_loglh!(c::ParticleCloud, loglh::Array{Float64,1})
    for (p,l) in zip(c.particles,loglh)
        p.loglh = l
    end
end

function update_logpost!(c::ParticleCloud, logpost::Array{Float64,1})
    for (p,l) in zip(c.particles,logpost)
        p.logpost = l
    end
end

function normalize_weights!(c::ParticleCloud)
    new_weights = get_weights(c)/sum(get_weights(c))
    for (p,w) in zip(c.particles, new_weights)
        update_weight!(p,w)
    end
end

function reset_weights!(c::ParticleCloud)
    map(p->update_weight!(p,1/length(c)),c.particles)
end

function update_mutation!(p::Particle, para::Array{Float64,1},
                          like::Float64, post::Float64, accept::Bool)
    p.value = para
    p.loglh = like
    p.logpost = post
    p.accept = accept
end

function update_val!(p::Particle, val::Array{Float64,1})
    p.value = val
end

function update_weight!(p::Particle, weight::Float64)
    p.weight = weight
end

function update_acceptance_rate!(c::ParticleCloud)
    temp = map(p -> p.accept, c.particles)
    arr = temp[1]
    for i in 2:length(temp)
        arr = vcat(arr,temp[i])
    end
    c.accept = mean(temp)
end

@inline function Base.length(c::ParticleCloud)
        return length(c.particles)
end

@inline function Base.length(p::Particle)
    return length(p.value)
end

function weighted_mean(c::ParticleCloud)
    return sum(map(*, c.particles, get_weights(c)))
end

function weighted_std(c::ParticleCloud)
    temp = map(p -> (p - weighted_mean(c)).^2, c.particles)
    return sqrt(sum(map(*, temp, get_weights(c))))
end

# Calculate the covariance of the particles
function Base.cov(c::ParticleCloud)
    temp = zeros(length(c.particles[1]),length(c))
    for (i,p) in enumerate(c.particles)
        temp[:,i] = p.value
    end
    return cov(temp,2)
end

for op in (:(Base.:+),
           :(Base.:-),
           :(Base.:*),
           :(Base.:/),
           :(Base.:^))

    @eval ($op)(p::Particle, q::Particle)   = ($op)(p.value,q.value)
    @eval ($op)(p::Particle, x::Integer)    = ($op)(p.value,x)
    @eval ($op)(p::Particle, x::Number)     = ($op)(p.value,x)
    @eval ($op)(x::Integer, p::Particle)    = ($op)(p.value,x)
    @eval ($op)(x::Number, p::Particle)     = ($op)(p.value,x)
    @eval ($op)(p::Particle, q::Array{Float64,1})   = ($op)(p.value,q)
    @eval ($op)(p::Array{Float64,1}, q::Particle)   = ($op)(p,q.value)
end
