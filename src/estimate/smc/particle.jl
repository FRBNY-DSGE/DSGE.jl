"""
```
Particle
```

The `Particle` type contains the values and weight of a given vector of parameters.

### Fields
- `weights::Vector{Float64}`: The weights of a given particle (from the perspective of importance sampling). The incremental weight that updates the weights at each step φ_n is given by p(Y|θ^i_{n-1})^(φ_n-φ_{n-1}).
- `keys::Vector{Symbol}`: The key of the parameter corresponding to a given value.
- `value::Vector{Float64}`: The matrix of particles (n_params by n_parts)
- `loglh::Float64`: The log-likelihood of the Particle
- `logpost::Float64`: The log-posterior of the Particle
- `old_loglh::Float64`: The log-likelihood of the Particle evaluated at the old data (only non-zero during time tempering)
- `old_logpost::Float64`: The log-posterior of the Particle evaluated at the old data (only non-zero during time tempering)
"""
type Particle
    weight::Float64
    keys::Vector{Symbol}
    value::Vector{Float64}
    loglh::Float64
    logpost::Float64
    old_loglh::Float64
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
- `particles::Vector{Particle}`: The vector of particles (which contain weight, keys, and value information)
- `tempering_schedule::Vector{Float64}`: The vector of ϕ_ns (tempering factors)
- `ESS::Vector{Float64}`: The vector of effective sample sizes (resample if ESS falls under the threshold)
- `stage_index::Int`: The current iteration index of the algorithm
- `n_Φ::Int`: The total number of stages of in the fixed tempering schedule
    (if the algorithm is run with an adaptive ϕ schedule then this is used to calibrate the ϕ_prop)
- `resamples::Int`: The number of times the particle population was resampled
- `c::Float64`: The mutation step size
- `accept::Float64`: The average acceptance rate of mutation steps
- `total_sampling_time::Float64`: The total amount of time that the smc algorith took to execute
"""
type ParticleCloud
    particles::Vector{Particle}
    tempering_schedule::Vector{Float64}
    ESS::Vector{Float64}
    stage_index::Int
    n_Φ::Int
    resamples::Int
    c::Float64
    accept::Float64
    total_sampling_time::Float64
end

# Easier constructor for ParticleCloud, which initializes the weights to be equal, and everything else
# in the Particle object etc. to be empty
function ParticleCloud(m::AbstractModel, n_parts::Int)
    return ParticleCloud([Particle(1/n_parts,[m.parameters[i].key for i in 1:length(m.parameters)],
                         zeros(length(m.parameters)),0.,0.,0.,false) for n in 1:n_parts],
                         zeros(1),zeros(1),1,0,0,0.,0., 0.)
end

function get_weights(c::ParticleCloud)
    return map(p -> p.weight,c.particles)
end

function get_vals(c::ParticleCloud)
    return hcat(map(p -> p.value, c.particles)...)
end

function get_loglh(c::ParticleCloud)
    return map(p -> p.loglh, c.particles)
end

function get_old_loglh(c::ParticleCloud)
    return map(p -> p.old_loglh, c.particles)
end

function get_logpost(c::ParticleCloud)
    return map(p -> p.logpost, c.particles)
end

function get_likeliest_particle_value(c::ParticleCloud)
    return c.particles[indmax(get_loglh(c))].value
end

function update_draws!(c::ParticleCloud, draws::Array{Float64,2})
    @assert size(draws) == (length(c.particles[1]),length(c))
    for (i,p) in enumerate(c.particles)
        update_val!(p, draws[:,i])
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
                          like::Float64, post::Float64, old_like::Float64, accept::Bool)
    p.value = para
    p.loglh = like
    p.logpost = post
    p.old_loglh = old_like
    p.accept = accept
end

function update_val!(p::Particle, val::Array{Float64,1})
    p.value = val
end

function update_weight!(p::Particle, weight::Float64)
    p.weight = weight
end

# For resetting a previously used cloud's settings for the purpose of time tempering
function reset_cloud_settings!(c::ParticleCloud)
    c.tempering_schedule = zeros(1)
    c.ESS = zeros(1)
    c.stage_index = 1
    c.n_Φ = 0
    c.resamples = 0
    c.c = 0
    c.accept = 0.
    c.total_sampling_time = 0.
    reset_weights!(c)
end

function update_acceptance_rate!(c::ParticleCloud)
    c.accept = mean(map(p -> p.accept, c.particles))
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
    return sqrt.(sum(map(*, temp, get_weights(c))))
end

# Calculate the covariance of the particles
function Base.cov(c::ParticleCloud)
    return cov(get_vals(c), 2)
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
