"""
```
one_draw(m::AbstractDSGEModel, data::Matrix{Float64}; use_chand_recursion::Bool = true,
         verbose::Symbol = :low)
```
Finds and returns one valid draw from parameter distribution, along with its log likelihood and log posterior.
"""
function one_draw(m::AbstractDSGEModel, data::Matrix{Float64};
                  use_chand_recursion::Bool = true, verbose::Symbol = :low)
    success    = false
    draw       = vec(rand(m.parameters, 1))
    draw_loglh = draw_logpost = 0.0

    while !success
        try
            update!(m, draw)
            draw_loglh = likelihood(m, data, catch_errors = true,
                                    use_chand_recursion = use_chand_recursion,
                                    verbose = verbose)
            draw_logpost = prior(m)
            if (draw_loglh == -Inf) | (draw_loglh === NaN)
                draw_loglh = draw_logpost = -Inf
            end
        catch err
            if isa(err, ParamBoundsError)
                draw_loglh = draw_logpost = -Inf
            elseif isa(err, PosDefException) || isa(err, SingularException) ||
                   isa(err, LinearAlgebra.LAPACKException) || isa(err, DSGE.SteadyStateConvergenceError)
                draw_loglh = draw_logpost = -Inf
            else
                throw(err)
            end
        end
        if any(isinf.(draw_loglh))
            draw = vec(rand(m.parameters, 1))
        else
            success = true
        end
    end
    return vector_reshape(draw, draw_loglh, draw_logpost)
end

"""
```
initial_draw!(m::AbstractDSGEModel, data::Matrix{Float64}, c::ParticleCloud)
initial_draw!(m::AbstractDSGEModel, data::Matrix{Float64}, c::Cloud)
```

Draw from a general starting distribution (set by default to be from the prior) to
initialize the SMC algorithm. Returns a tuple (logpost, loglh) and modifies the
particle objects in the particle cloud in place.
"""
function initial_draw!(m::AbstractDSGEModel, data::Matrix{Float64},
                       c::Union{Cloud, ParticleCloud};
                       parallel::Bool = false, use_chand_recursion::Bool = true,
                       verbose::Symbol = :low)
    n_parts = length(c)

    # ================== Define closure on one_draw function ==================
    sendto(workers(), m = m)
    sendto(workers(), data = data)
    sendto(workers(), verbose = verbose)
    sendto(workers(), use_chand_recursion = use_chand_recursion)

    one_draw_closure() = one_draw(m, data; use_chand_recursion = use_chand_recursion,
                                  verbose = verbose)
    @everywhere one_draw_closure() = one_draw(m, data;
                                              use_chand_recursion = use_chand_recursion,
                                              verbose = verbose)
    # =========================================================================

    # For each particle, finds valid parameter draw and returns likelihood & posterior
    draws, loglh, logpost = if parallel
        @sync @distributed (vector_reduce) for i in 1:n_parts
            one_draw_closure()
        end
    else
        vector_reduce([one_draw_closure() for i in 1:n_parts]...)
    end

    update_draws!(c, draws)
    update_loglh!(c, vec(loglh))
    update_logpost!(c, vec(logpost))
    update_old_loglh!(c, zeros(n_parts))
    # Need to call `set_weights` as opposed to `update_weights`
    # since update_weights will multiply and 0*anything = 0
    set_weights!(c, ones(n_parts))
end

"""
```
function draw_likelihood(m, data, draw_vec; verbose::Symbol = :low)
```
Computes likelihood of a particular parameter draw; returns loglh and logpost.
"""
function draw_likelihood(m::AbstractDSGEModel, data::Matrix{Float64},
                         draw_vec::Vector{Float64}; verbose::Symbol = :low)
    update!(m, draw_vec)
    loglh   = likelihood(m, data, verbose = verbose)
    logpost = prior(m)
    return scalar_reshape(loglh, logpost)
end

"""
```
initialize_likelihoods!(m::AbstractDSGEModel, data::Matrix{Float64},
                        c::Union{Cloud, ParticleCloud};
                        parallel::Bool = false, verbose::Symbol = :low)
```
This function is made for transfering the log-likelihood values saved in the
Cloud from a previous estimation to each particle's respective old_loglh
field, and for evaluating/saving the likelihood and posterior at the new data, which
here is just the argument, data.
"""
function initialize_likelihoods!(m::AbstractDSGEModel, data::Matrix{Float64},
                                 c::Union{Cloud, ParticleCloud};
                                 parallel::Bool = false, verbose::Symbol = :low)
    n_parts = length(c)
    draws = (typeof(c) <: Cloud) ? get_vals(c; transpose = false) : Matrix{Float64}(get_vals(c)')

    # Retire log-likelihood values from the old estimation to the field old_loglh
    update_old_loglh!(c, get_loglh(c))

    # ============== Define closure on draw_likelihood function ===============
    sendto(workers(), m = m)
    sendto(workers(), data = data)
    sendto(workers(), verbose = verbose)

    draw_likelihood_closure(draw::Vector{Float64}) = draw_likelihood(m, data, draw;
                                                                 verbose = verbose)
    @everywhere draw_likelihood_closure(draw::Vector{Float64}) = draw_likelihood(m, data,
                                                                 draw; verbose = verbose)
    # =========================================================================

    # TODO: handle when the likelihood with new data cannot be evaluated (returns -Inf),
    # even if the likelihood was not -Inf prior to incorporating new data
    loglh, logpost = if parallel
        @sync @distributed (scalar_reduce) for i in 1:n_parts
            draw_likelihood_closure(draws[i, :])
        end
    else
        scalar_reduce([draw_likelihood_closure(draws[i, :]) for i in 1:n_parts]...)
    end
    update_loglh!(c, loglh)
    update_logpost!(c, logpost)
end

"""
```
function initialize_cloud_settings!(m::AbstractDSGEModel, cloud::ParticleCloud;
                                    tempered_update::Bool = false)
```
Initializes stage index, number of Φ stages, c, resamples, acceptance, and sampling time.
"""
function initialize_cloud_settings!(m::AbstractDSGEModel,
                                    cloud::Union{ParticleCloud,Cloud};
                                    tempered_update::Bool = false)
    if tempered_update
        cloud.ESS = [cloud.ESS[end]]
    else
        cloud.ESS[1] = get_setting(m, :n_particles)
    end
    cloud.stage_index = 1
    cloud.n_Φ         = get_setting(m, :n_Φ)
    cloud.resamples   = 0
    cloud.c           = get_setting(m, :step_size_smc)
    cloud.accept      = get_setting(m, :target_accept)
    cloud.total_sampling_time = 0.
    cloud.tempering_schedule  = zeros(1)
end
