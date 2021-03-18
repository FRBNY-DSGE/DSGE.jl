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
- `accept::Float64`: Acceptance rate of Metropolis-Hastings steps in the mutation of the particle in the previous period
"""
mutable struct Particle
    weight::Float64
    keys::Vector{Symbol}
    value::Vector{Float64}
    loglh::Float64
    logpost::Float64
    old_loglh::Float64
    accept::Float64
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

# Implementing a sensible showing option for clouds within DataFrames
function Base.show(io::IO, cloud::ParticleCloud)
    @printf io "Particle Cloud: N_parts=%i, N_stages=%i" length(cloud) cloud.stage_index
end

"""
```
function ParticleCloud(m::AbstractDSGEModel, n_parts::Int)
function ParticleCloud(m::AbstractDSGEVARModel, n_parts::Int)
```
Easier constructor for ParticleCloud, which initializes the weights to be equal, and everything else in the Particle object etc. to be empty.
"""
function ParticleCloud(m::AbstractDSGEModel, n_parts::Int)
    return ParticleCloud([Particle(1/n_parts,[m.parameters[i].key for i in 1:length(m.parameters)],
                         zeros(length(m.parameters)),0.,0.,0.,false) for n in 1:n_parts],
                         zeros(1),zeros(1),1,0,0,0.,0.25, 0.)
end

function ParticleCloud(m::AbstractDSGEVARModel, n_parts::Int)
    return ParticleCloud(m.dsge, n_parts)
end


"""
```
Cloud
```

The `Cloud` type contains all of the relevant information for a given cloud of particles
in the SMC algorithm. Information for a single iteration is stored at any given time (and thus
the final output will be the final cloud of particles, of which only the particle values will be saved).

### Fields
- `particles::Matrix{Float64}`: The vector of particles (which contain weight, keys, and value information)
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
mutable struct Cloud
    particles::Matrix{Float64}
    tempering_schedule::Vector{Float64}
    ESS::Vector{Float64}
    stage_index::Int
    n_Φ::Int
    resamples::Int
    c::Float64
    accept::Float64
    total_sampling_time::Float64
end

"""
```
function Cloud(m::AbstractDSGEModel, n_parts::Int)
```
Easier constructor for Cloud, which initializes the weights to be equal, and everything else in the Particle object etc. to be empty.
"""
function Cloud(m::AbstractDSGEModel, n_parts::Int)
    return Cloud(Matrix{Float64}(undef, n_parts, n_parameters(m) + 6),
                 zeros(1),zeros(1),1,0,0,0.,0.25, 0.)
end

"""
```
function Cloud(c::ParticleCloud)
function Cloud(c::Cloud)
```
Returns a Cloud type.
"""
function Cloud(cloud::ParticleCloud)
    return Cloud(hcat(Matrix{Float64}(get_vals(cloud)'), get_loglh(cloud), get_logprior(cloud),
                      get_old_loglh(cloud), get_accept(cloud), get_weights(cloud)),
                 cloud.tempering_schedule, cloud.ESS, cloud.stage_index, cloud.n_Φ,
                 cloud.resamples, cloud.c, cloud.accept, cloud.total_sampling_time)
end
function Cloud(cloud::Cloud)
    return cloud
end


"""
```
function ParticleCloud(cloud::Cloud, para_symbols::Vector{Symbol})
```
Returns a ParticleCloud in the place of a Cloud.
"""
function ParticleCloud(cloud::Union{SMC.Cloud, DSGE.Cloud}, para_symbols::Vector{Symbol})
    N = size(cloud.particles, 2)
    return ParticleCloud([Particle(cloud.particles[k, ind_weight(N)], para_symbols,
                                   cloud.particles[k, 1:ind_para_end(N)],
                                   cloud.particles[k, ind_loglh(N)],
                                   cloud.particles[k, ind_logpost(N)], # use logpost here b/c column index is typically retrieved via ind_logpost for ParticleCloud
                                   cloud.particles[k, ind_old_loglh(N)],
                                   cloud.particles[k, ind_accept(N)]) for k in 1:length(cloud)],
                         cloud.tempering_schedule, cloud.ESS, cloud.stage_index, cloud.n_Φ,
                         cloud.resamples, cloud.c, cloud.accept, cloud.total_sampling_time)
end

"""
Find correct indices for accessing columns of cloud array.
"""
ind_para_end(N::Int)  = N-5
ind_loglh(N::Int)     = N-4
ind_logpost(N::Int)   = N-3
ind_logprior(N::Int)  = N-3 # TODO: Fix logprior/logpost shenanigans. Basically, logpost really should be logprior b/c the logprior is what's stored.
ind_old_loglh(N::Int) = N-2
ind_accept(N::Int)    = N-1
ind_weight(N::Int)    = N

"""
```
function get_weights(c::ParticleCloud)
function get_weights(c::Cloud)
function get_weights(c::Matrix{Float64})
```
Returns Vector{Float64}(n_parts) of weights of particles in cloud.
"""
function get_weights(c::ParticleCloud)
    return map(p -> p.weight, c.particles)
end
@inline function get_weights(c::Matrix{Float64})
    return c[:, ind_weight(size(c,2))]
end
@inline function get_weights(c::Union{DSGE.Cloud, SMC.Cloud})
    return c.particles[:, ind_weight(size(c.particles,2))]
end

"""
```
function get_vals(c::ParticleCloud)
function get_vals(c::Matrix{Float64})
function get_vals(c::Cloud)
```
Returns Matrix{Float64}(n_params, n_parts) of parameter values in particle cloud.
"""
function get_vals(c::ParticleCloud)
    return hcat(map(p -> p.value, c.particles)...)
end
@inline function get_vals(c::Matrix{Float64}; transpose::Bool = true)
    return transpose ? Matrix{Float64}(c[:, 1:ind_para_end(size(c, 2))]') :
                                       c[:, 1:ind_para_end(size(c, 2))]
end
@inline function get_vals(c::Cloud; transpose::Bool = true)
    return transpose ? Matrix{Float64}(c.particles[:, 1:ind_para_end(size(c.particles,2))]') :
                                       c.particles[:, 1:ind_para_end(size(c.particles,2))]
end

"""
```
function get_loglh(c::ParticleCloud)
function get_loglh(c::Matrix{Float64})
function get_loglh(c::Cloud)
```
Returns Vector{Float64}(n_parts) of log-likelihood of particles in cloud.
"""
function get_loglh(c::ParticleCloud)
    return map(p -> p.loglh, c.particles)
end
@inline function get_loglh(c::Matrix{Float64})
    return c[:, ind_loglh(size(c,2))]
end
@inline function get_loglh(c::Cloud)
    return c.particles[:, ind_loglh(size(c.particles,2))]
end

"""
```
function cloud_isempty(c::ParticleCloud)
function cloud_isempty(c::Matrix{Float64})
function cloud_isempty(c::Cloud)
```
Check if cloud has no particles.
"""
function cloud_isempty(c::ParticleCloud)
    return isempty(c.particles)
end
function cloud_isempty(c::Matrix{Float64})
    return isempty(c)
end
function cloud_isempty(c::Cloud)
    println("TODO: 'ISEMPTY' IS INCOMPATIBLE WITH INITIALIZATION")
    return isempty(c.particles)
end


"""
```
function get_old_loglh(c::ParticleCloud)
function get_old_loglh(c::Matrix{Float64})
function get_old_loglh(c::Cloud)
```
Returns Vector{Float64}(n_parts) of old log-likelihood of particles in cloud.
"""
function get_old_loglh(c::ParticleCloud)
    return map(p -> p.old_loglh, c.particles)
end
@inline function get_old_loglh(c::Matrix{Float64})
    return c[:, ind_old_loglh(size(c,2))]
end
@inline function get_old_loglh(c::Cloud)
    return c.particles[:, ind_old_loglh(size(c.particles,2))]
end

"""
```
function get_logpost(c::ParticleCloud)
function get_logpost(c::Matrix{Float64})
function get_logpost(c::Cloud)
```
Returns Vector{Float64}(n_parts) of log-posterior of particles in cloud.
"""
function get_logpost(c::ParticleCloud)
    return map(p -> p.logpost + p.loglh, c.particles)
end
@inline function get_logpost(c::Matrix{Float64})
    return c[:, ind_loglh(size(c,2))] .+ c[:, ind_logprior(size(c,2))]
end
@inline function get_logpost(c::Cloud)
    return c.particles[:, ind_loglh(size(c.particles,2))] .+
        c.particles[:, ind_logprior(size(c.particles,2))]
end

"""
```
function get_logprior(c::ParticleCloud)
function get_logprior(c::Matrix{Float64})
function get_logprior(c::Cloud)
```
Returns Vector{Float64}(n_parts) of log-prior of particles in cloud.
"""
function get_logprior(c::ParticleCloud)
    #TODO: Fix logpost/logprior confusion
    return map(p -> p.logpost, c.particles)
end
@inline function get_logprior(c::Matrix{Float64})
    return c[:, ind_logprior(size(c,2))]
end
@inline function get_logprior(c::Cloud)
    return c.particles[:, ind_logprior(size(c.particles,2))]
end

"""
```
function get_accept(c::ParticleCloud)
function get_accept(c::Matrix{Float64})
function get_accept(c::Cloud)
```
Returns Vector{Float64}(n_parts) of old log-likelihood of particles in cloud.
"""
function get_accept(c::ParticleCloud)
    return map(p -> p.accept, c.particles)
end
@inline function get_accept(c::Matrix{Float64})
    return c[:, ind_accept(size(c,2))]
end
@inline function get_accept(c::Cloud)
    return c.particles[:, ind_accept(size(c.particles,2))]
end

"""
```
function get_likeliest_particle_value(c::ParticleCloud)
function get_likeliest_particle_value(c::Matrix{Float64})
function get_likeliest_particle_value(c::Cloud)
```
Return parameter vector of particle with highest log-likelihood.
"""
function get_likeliest_particle_value(c::ParticleCloud)
    return c.particles[indmax(get_loglh(c))].value
end
@inline function get_likeliest_particle_value(c::Matrix{Float64})
    return c[indmax(get_loglh(c)), 1:ind_para_end(size(c,2))]
end
@inline function get_likeliest_particle_value(c::Cloud)
    return c.particles[indmax(get_loglh(c)), 1:ind_para_end(size(c.particles,2))]
end

"""
```
function update_draws!(c::ParticleCloud, draws::Matrix{Float64})
function update_draws!(c::ParticleCloud, draws::Array{Particle,1})
@inline function update_draws!(c::Matrix{Float64}, draws::Matrix{Float64})
```
Update parameter draws in cloud.
"""
function update_draws!(c::ParticleCloud, draws::Matrix{Float64})
    @assert size(draws) == (length(c.particles[1]), length(c))
    for (i,p) in enumerate(c.particles)
        update_val!(p, draws[:,i])
    end
end
function update_draws!(c::ParticleCloud, draws::Array{Particle,1})
    for (i,d) in enumerate(draws)
        c.particles[i] = d
    end
end
@inline function update_draws!(c::Cloud, draws::Matrix{Float64})
    I, J     = size(draws)
    n_parts  = length(c)
    n_params = ind_para_end(size(c.particles, 2))
    if (I, J) == (n_parts, n_params)
        for i = 1:I, j=1:J
            c.particles[i, j] = draws[i, j]
        end
    elseif (I, J) == (n_params, n_parts)
        for i = 1:I, j=1:J
            c.particles[j, i] = draws[i, j]
        end
    else
        throw(error("update_draws!(c::Cloud, draws::Matrix): Draws are incorrectly sized!"))
    end
end

"""
```
function update_weights!(c::ParticleCloud, incweight::Vector{Float64})
@inline function update_weights!(c::Matrix{Float64}, incweight::Vector{Float64})
@inline function update_weights!(c::Cloud, weights::Vector{Float64})
```
Update weights in cloud.
"""
function update_weights!(c::ParticleCloud, incweight::Vector{Float64})
    new_weights = get_weights(c) .* incweight
    for (p,w) in zip(c.particles, new_weights)
        update_weight!(p, w)
    end
end
@inline function update_weights!(c::Matrix{Float64}, incweight::Vector{Float64})
    @assert size(c, 1) == length(incweight) "Dimensional mismatch in inc. weights"
    N = ind_weight(size(c,2))
    for i=1:length(incweight)
        c[i, N] *= incweight[i]
    end
end
@inline function update_weights!(c::Cloud, weights::Vector{Float64})
    update_weights!(c.particles, weights)
end

"""
```
@inline function set_weights!(c::Cloud, weights::Vector{Float64})
```
Set weights to specific values. Contrast to update_weights, which multiplies
existing weights by provided incremental weights.
"""
function set_weights!(c::ParticleCloud, weights::Vector{Float64})
    for (p,w) in zip(c.particles, weights)
        update_weight!(p,w)
    end
end
@inline function set_weights!(c::Cloud, weights::Vector{Float64})
    @assert length(c) == length(weights) "Dimensional mismatch in set_weights"
    N = ind_weight(size(c.particles,2))
    for i=1:length(c)
        c.particles[i, N] = weights[i]
    end
end


"""
```
function update_loglh!(c::ParticleCloud, incweight::Vector{Float64})
@inline function update_loglh!(c::Matrix{Float64}, incweight::Vector{Float64})
@inline function update_loglh!(c::Cloud, loglh::Vector{Float64})
```
Update log-likelihood in cloud.
"""
function update_loglh!(c::ParticleCloud, loglh::Vector{Float64})
    for (p,l) in zip(c.particles,loglh)
        p.loglh = l
    end
end
@inline function update_loglh!(c::Matrix{Float64}, loglh::Vector{Float64})
    @assert size(c,1) == length(loglh) "Dimensional mismatch"
    N = ind_loglh(size(c,2))
    for i=1:length(loglh)
        c[i, N] = loglh[i]
    end
end
@inline function update_loglh!(c::Cloud, loglh::Vector{Float64})
    update_loglh!(c.particles, loglh)
end

"""
```
function update_logpost!(c::ParticleCloud, incweight::Vector{Float64})
@inline function update_logpost!(c::Matrix{Float64}, incweight::Vector{Float64})
@inline function update_logpost!(c::Cloud, logpost::Vector{Float64})
```
Update log-posterior in cloud.
"""
function update_logpost!(c::ParticleCloud, logpost::Vector{Float64})
    for (p,l) in zip(c.particles,logpost)
        p.logpost = l
    end
end
@inline function update_logpost!(c::Matrix{Float64}, logpost::Vector{Float64})
    @assert size(c, 1) == length(logpost) "Dimensional mismatch"
    N = ind_logpost(size(c,2))
    for i=1:length(logpost)
        c[i, N] = logpost[i]
    end
end
@inline function update_logpost!(c::Cloud, logpost::Vector{Float64})
    update_logpost!(c.particles, logpost)
end


"""
```
function update_old_loglh!(c::ParticleCloud, incweight::Vector{Float64})
@inline function update_old_loglh!(c::Matrix{Float64}, incweight::Vector{Float64})
@inline function update_old_loglh!(c::Cloud, old_loglh::Vector{Float64})
```
Update log-likelihood in cloud.
"""
function update_old_loglh!(c::ParticleCloud, old_loglh::Vector{Float64})
    for (p,l) in zip(c.particles,old_loglh)
        p.old_loglh = l
    end
end
@inline function update_old_loglh!(c::Matrix{Float64}, old_loglh::Vector{Float64})
    @assert size(c, 1) == length(old_loglh) "Dimensional mismatch"
    N = ind_old_loglh(size(c,2))
    for i=1:length(old_loglh)
        c[i, N] = old_loglh[i]
    end
end
@inline function update_old_loglh!(c::Cloud, old_loglh::Vector{Float64})
    update_old_loglh!(c.particles, old_loglh)
end

"""
```
function normalize_weights!(c::ParticleCloud)
function normalize_weights!(c::Matrix{Float64})
@inline function normalize_weights!(c::Cloud)
```
Normalize weights in cloud.
"""
function normalize_weights!(c::ParticleCloud)
    new_weights = get_weights(c) / sum(get_weights(c))
    for (p,w) in zip(c.particles, new_weights)
        update_weight!(p,w)
    end
end
@inline function normalize_weights!(c::Matrix{Float64})
    sum_weights = sum(get_weights(c))
    c[:, ind_weight(size(c,2))] /= sum_weights
end
@inline function normalize_weights!(c::Cloud)
    normalize_weights!(c.particles)
end

"""
```
function reset_weights!(c::ParticleCloud)
function reset_weights!(c::Matrix{Float64})
function reset_weights!(c::Cloud)
```
Uniformly reset weights of all particles to 1/n_parts.
"""
function reset_weights!(c::ParticleCloud)
    map(p->update_weight!(p,1/length(c)), c.particles)
end
@inline function reset_weights!(c::Matrix{Float64})
    n_parts = size(c, 1)
    c[:, ind_weight(size(c,2))] .= 1.0 / n_parts
end
function reset_weights!(c::Cloud)
    reset_weights!(c.particles)
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

"""
```
function update_mutation!(p::Particle, para::Vector{Float64},
                          like::Float64, post::Float64, old_like::Float64, accept::Float64)
function update_mutation!(p::Vector{Float64}, para::Vector{Float64},
                          like::Float64, post::Float64, old_like::Float64, accept::Float64)
```
Update a particle's parameter vector, log-likelihood, log-posteriod, old log-likelihood, and acceptance rate at the end of mutation.
"""
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
    N = length(p)
    p[1:ind_para_end(N)] = para
    p[ind_loglh(N)]      = like
    p[ind_logpost(N)]    = post
    p[ind_old_loglh(N)]  = old_like
    p[ind_accept(N)]     = accept
end

"""
```
function update_cloud!(cloud::ParticleCloud, new_particles::Matrix{Float64})
function update_cloud!(cloud::Cloud, new_particles::Matrix{Float64})
```
Updates cloud values with those of new particles in particle array.
"""
function update_cloud!(cloud::ParticleCloud, new_particles::Matrix{Float64})
    M, N = size(new_particles)
    if M == length(cloud.particles)
        for k=1:length(cloud.particles)
            cloud.particles[k].value     = new_particles[k, 1:ind_para_end(N)]
            cloud.particles[k].loglh     = new_particles[k, ind_loglh(N)]
            cloud.particles[k].logpost   = new_particles[k, ind_logpost(N)]
            cloud.particles[k].old_loglh = new_particles[k, ind_old_loglh(N)]
            cloud.particles[k].accept    = new_particles[k, ind_accept(N)]
        end
    elseif N == length(cloud.particles)
        for k=1:length(cloud.particles)
            cloud.particles[k].value     = new_particles[1:ind_para_end(M), k]
            cloud.particles[k].loglh     = new_particles[ind_loglh(M), k]
            cloud.particles[k].logpost   = new_particles[ind_logpost(M), k]
            cloud.particles[k].old_loglh = new_particles[ind_old_loglh(M), k]
            cloud.particles[k].accept    = new_particles[ind_accept(M), k]
        end
    else
        throw(error("update_cloud!(c::ParticleCloud, draws): draws are incorrect size!"))
    end
end
function update_cloud!(cloud::Cloud, new_particles::Matrix{Float64})
    I, J = size(new_particles)
    if I == length(cloud)
        cloud.particles = new_particles
    elseif J == length(cloud)
        for k = 1:length(cloud)
            cloud.particles[k,:] = new_particles[:, k]
        end
    else
        throw(error("update_cloud!(c::Cloud, draws): draws are incorrect size!"))
    end
end


"""
```
function get_cloud_array(c::ParticleCloud)
```
Returns a particle array in the place of a ParticleCloud.
"""
function get_cloud_array(c::ParticleCloud)
    return hcat(Matrix{Float64}(get_vals(c)'), get_loglh(c), get_logprior(c),
                get_old_loglh(c), get_accept(c), get_weights(c))
end

"""
```
function update_val!(p::Particle, val::Vector{Float64})
@inline function update_val!(p::Vector{Float64}, val::Vector{Float64})
```
Update parameter vector of particle.
"""
function update_val!(p::Particle, val::Vector{Float64})
    p.value = val
end
@inline function update_val!(p::Vector{Float64}, val::Vector{Float64})
    @assert ind_para_end(length(p)) == length(val) "Parameter vector length is wrong!"
    p[1:length(val)] = val
end

"""
```
function update_weight!(p::Particle, weight::Vector{Float64})
@inline function update_weight!(p::Vector{Float64}, weight::Vector{Float64})
```
Update weight of particle.
"""
function update_weight!(p::Particle, weight::Float64)
    p.weight = weight
end
@inline function update_weight!(p::Vector{Float64}, weight::Float64)
    p[ind_weight(length(p))] = weight
end

"""
```
function update_acceptance_rate!(c::ParticleCloud)
function update_acceptance_rate!(c::Cloud)
```
Update cloud's acceptance rate with the mean of its particle acceptance rates.
"""
function update_acceptance_rate!(c::ParticleCloud)
    c.accept = mean(map(p -> p.accept, c.particles))
end
function update_acceptance_rate!(c::Cloud)
    c.accept = mean(get_accept(c))
end

@inline function Base.length(c::ParticleCloud)
        return length(c.particles)
end
@inline function Base.length(c::Cloud)
        return size(c.particles, 1)
end
@inline function Base.length(p::Particle)
    return length(p.value)
end

"""
```
function weighted_mean(c::ParticleCloud)
function weighted_mean(c::Cloud)
function weighted_mean(c::Matrix{Float64})
```
Compute weighted mean of particle cloud.
"""
function weighted_mean(c::ParticleCloud)
    return get_vals(c) * get_weights(c)
end
@inline function weighted_mean(c::Matrix{Float64})
    return get_vals(c) * get_weights(c)
end
@inline function weighted_mean(c::Cloud)
    return weighted_mean(c.particles)
end

"""
```
function weighted_quantile(c::ParticleCloud, i::Int64)
@inline function weighted_quantile(c::Matrix{Float64}, i::Int64)
function weighted_quantile(c::Cloud, i::Int64)
```
Compute weighted quantiles of particle cloud for input parameter, indexed by i.
"""
function weighted_quantile(c::ParticleCloud, i::Int64)
    lb = quantile(get_vals(c)[i, :], Weights(get_weights(c)), .05)
    ub = quantile(get_vals(c)[i, :], Weights(get_weights(c)), .95)
    return lb, ub
end
@inline function weighted_quantile(c::Matrix{Float64}, i::Int64)
    @assert i <= ind_end_para(size(c,2)) "Parameter index invalid."
    lb = quantile(c[:, i], Weights(get_weights(c)), .05)
    ub = quantile(c[:, i], Weights(get_weights(c)), .95)
    return lb, ub
end
function weighted_quantile(c::Cloud, i::Int64)
    return weighted_quantile(c.particles, i)
end

"""
```
function weighted_std(c::ParticleCloud)
function weighted_std(c::Cloud)
@inline function weighted_std(c::Matrix{Float64})
```
Compute weighted standard deviation of particle cloud.
"""
function weighted_std(c::ParticleCloud)
    return sqrt.(diag(weighted_cov(c)))
end
@inline function weighted_std(c::Matrix{Float64})
    return sqrt.(diag(weighted_cov(c)))
end
function weighted_std(c::Cloud)
    return weighted_std(c.particles)
end

"""
```
function weighted_cov(c::ParticleCloud)
function weighted_cov(c::Cloud)
@inline function weighted_cov(c::Matrix{Float64})
```
Compute weighted covariance of particle cloud.
"""
function weighted_cov(c::ParticleCloud)
    return cov(Matrix(get_vals(c)'), Weights(get_weights(c)), corrected = false)
end
@inline function weighted_cov(c::Matrix{Float64})
    return cov(get_vals(c; transpose = false), Weights(get_weights(c)), corrected = false)
end
@inline function weighted_cov(c::Cloud)
    return weighted_cov(c.particles)
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
    @eval ($op)(p::Particle, q::Vector{Float64})   = @. ($op)(p.value,q)
    @eval ($op)(p::Vector{Float64}, q::Particle)   = @. ($op)(p,q.value)
end
