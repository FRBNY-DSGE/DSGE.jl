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
- `accept::Bool`: Whether or not the mutation of the particle was accepted in the previous period
"""
mutable struct Particle
    weight::Float64
    keys::Vector{Symbol}
    value::Vector{Float64}
    loglh::Float64
    logpost::Float64
    old_loglh::Float64
    accept::Float64 #Bool # TO MATCH FORTRAN
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
mutable struct ParticleCloud
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
                         zeros(1),zeros(1),1,0,0,0.,0.25, 0.)
end

function get_weights(c::ParticleCloud)
    return map(p -> p.weight, c.particles)
end
@inline function get_weights(c::Matrix{Float64})
    return c[:,end]
end

function get_vals(c::ParticleCloud)
    return hcat(map(p -> p.value, c.particles)...)
end
@inline function get_vals(c::Matrix{Float64})
    return c[:, 1:end-5]
end

function get_loglh(c::ParticleCloud)
    return map(p -> p.loglh, c.particles)
end
@inline function get_loglh(c::Matrix{Float64})
    return c[:, 1:end-4]
end

function cloud_isempty(c::ParticleCloud)
    return isempty(c.particles)
end

function get_old_loglh(c::ParticleCloud)
    return map(p -> p.old_loglh, c.particles)
end
@inline function get_old_loglh(c::Matrix{Float64})
    return c[:,end-2]
end

function get_logpost(c::ParticleCloud)
    return map(p -> p.logpost + p.loglh, c.particles)
end
@inline function get_logpost(c::Matrix{Float64})
    return c[:,end-3] .+ c[:,end-4]
end

#TODO: Fix logpost/logprior confusion
function get_logpprior(c::ParticleCloud)
    return map(p -> p.logpost, c.particles)
end
@inline function get_logprior(c::Matrix{Float64})
    return c[:,end-3]
end

function get_likeliest_particle_value(c::ParticleCloud)
    return c.particles[indmax(get_loglh(c))].value
end
@inline function get_likeliest_particle_value(c::Matrix{Float64})
    return c[indmax(get_loglh(c)), 1:end-4]
end

function update_draws!(c::ParticleCloud, draws::Matrix{Float64})
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
@inline function update_draws!(c::Matrix{Float64}, draws::Matrix{Float64})
    @assert size(draws) == (size(c,1), size(c,2)-5) "Draws are incorrectly sized"
    I, J = size(draws)
    for i = 1:I, j=1:J
        c[i, j] = draws[i,j]
    end
end

function update_weights!(c::ParticleCloud, incweight::Vector{Float64})
    new_weights = get_weights(c).*incweight
    for (p,w) in zip(c.particles, new_weights)
        update_weight!(p, w)
    end
end
@inline function update_weights!(c::Matrix{Float64}, incweight::Vector{Float64})
    @assert size(c, 1) == length(incweight) "Dimensional mismatch in inc. weights"
    for i=1:length(incweight)
        c[i,end] *= incweight[i]
    end
end

function update_loglh!(c::ParticleCloud, loglh::Vector{Float64})
    for (p,l) in zip(c.particles,loglh)
        p.loglh = l
    end
end
@inline function update_loglh!(c::Matrix{Float64}, loglh::Vector{Float64})
    @assert size(c, 1) == length(loglh) "Dimensional mismatch"
    for i=1:length(loglh)
        c[i, end-4] = loglh[i]
    end
end

function update_logpost!(c::ParticleCloud, logpost::Vector{Float64})
    for (p,l) in zip(c.particles,logpost)
        p.logpost = l
    end
end
@inline function update_logpost!(c::Matrix{Float64}, logpost::Vector{Float64})
    @assert size(c, 1) == length(logpost) "Dimensional mismatch"
    for i=1:length(logpost)
        c[i, end-3] = logpost[i]
    end
end

function normalize_weights!(c::ParticleCloud)
    new_weights = get_weights(c) / sum(get_weights(c))
    for (p,w) in zip(c.particles, new_weights)
        update_weight!(p,w)
    end
end
@inline function normalize_weights!(c::Matrix{Float64})
    sum_weights = sum(get_weights(c))
    c[:, end] /= sum_weights
end

function reset_weights!(c::ParticleCloud)
    map(p->update_weight!(p,1/length(c)),c.particles)
end
@inline function reset_weights!(c::Matrix{Float64})
    n_parts = size(c, 1)
    c[:, end] .= 1.0 / n_parts
end

"""
```
function update_mutation!(p::Particle, para::Vector{Float64},
                          like::Float64, post::Float64, old_like::Float64, accept::Bool)
```
Deprecated syntax. Accept should be a Float64, not a Bool.
"""
function update_mutation!(p::Particle, para::Vector{Float64},
                          like::Float64, post::Float64, old_like::Float64, accept::Bool)
    p.value = para
    p.loglh = like
    p.logpost = post
    p.old_loglh = old_like
    p.accept = accept
end
function update_mutation!(p::Particle, para::Vector{Float64},
                          like::Float64, post::Float64, old_like::Float64, accept::Float64)
    p.value = para
    p.loglh = like
    p.logpost = post
    p.old_loglh = old_like
    p.accept = accept
end
@inline function update_mutation!(p::Vector{Float64}, para::Vector{Float64}, like::Float64,
                                  post::Float64, old_like::Float64, accept::Float64)
    p[1:end-5] = para
    p[end-4]   = like
    p[end-3]   = post
    p[end-2]   = old_like
    p[end-1]   = accept
end

function update_val!(p::Particle, val::Vector{Float64})
    p.value = val
end
@inline function update_val!(p::Vector{Float64}, val::Vector{Float64})
    @assert length(p)-5 == length(val) "Parameter vector length is wrong!"
    p[1:length(val)] = val
end

function update_weight!(p::Particle, weight::Float64)
    p.weight = weight
end
@inline function update_weight!(p::Vector{Float64}, weight::Float64)
    p[end] = weight
end

function update_acceptance_rate!(c::ParticleCloud)
    c.accept = mean(map(p -> p.accept, c.particles))
end
@inline function mean_acceptance_rate(c::Matrix{Float64})
    return mean(c[:,end-1])
end

# TODO: Check this never gets called on new Particle type
@inline function Base.length(c::ParticleCloud)
        return length(c.particles)
end
@inline function Base.length(p::Particle)
    return length(p.value)
end

function weighted_mean(c::ParticleCloud)
    return dropdims(mean(get_vals(c), Weights(get_weights(c)), dims = 2), dims = 2)
end
@inline function weighted_mean(c::Matrix{Float64})
    return dropdims(mean(get_vals(c), Weights(get_weights(c)), dims = 2), dims = 2)
end

function weighted_quantile(c::ParticleCloud, i::Int64)
    lb = quantile(get_vals(c)[i, :], Weights(get_weights(c)), .05)
    ub = quantile(get_vals(c)[i, :], Weights(get_weights(c)), .95)
    return lb, ub
end
@inline function weighted_quantile(c::Matrix{Float64}, i::Int64)
    lb = quantile(c[:, i], Weights(c[:,end]), .05)
    ub = quantile(c[:, i], Weights(c[:,end]), .95)
    return lb, ub
end

function weighted_std(c::ParticleCloud)
    return sqrt.(diag(weighted_cov(c)))
end
@inline function weighted_std(c::Matrix{Float64})
    return sqrt.(diag(weighted_cov(c)))
end

function weighted_cov(c::ParticleCloud)
    return cov(Matrix(get_vals(c)'), Weights(get_weights(c)), corrected = false)
end
@inline function weighted_cov(c::Matrix{Float64})
    return cov(Matrix(get_vals(c)), Weights(get_weights(c)), corrected = false)
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
    @eval ($op)(p::Particle, q::Vector{Float64})   = ($op)(p.value,q)
    @eval ($op)(p::Vector{Float64}, q::Particle)   = ($op)(p,q.value)
end
