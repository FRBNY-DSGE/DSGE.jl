"""
```
initial_draw!(m::AbstractModel, data::Matrix{Float64}, c::ParticleCloud)
```

Draw from a general starting distribution (set by default to be from the prior) to initialize the SMC algorithm.
Returns a tuple (logpost, loglh) and modifies the particle objects in the particle cloud in place.

"""
function initial_draw!(m::AbstractModel, data::AbstractMatrix, c::ParticleCloud;
                       parallel::Bool = false, use_chand_recursion::Bool = false, verbose::Symbol = :low)
    n_parts = length(c)
    loglh = zeros(n_parts)
    logpost = zeros(n_parts)
    if parallel
        draws, loglh, logpost = @distributed (vector_reduce) for i in 1:n_parts
            draw = vec(rand(m.parameters, 1))
            draw_loglh = 0.
            draw_logpost = 0.
            success = false
            while !success
                try
                    update!(m, draw)
                    draw_loglh   = likelihood(m, data, catch_errors = true, use_chand_recursion=use_chand_recursion, verbose = verbose)
                    draw_logpost = prior(m)
                    if (draw_loglh == -Inf) | (draw_loglh===NaN)
                        draw_logpost = -Inf
                        draw_loglh = -Inf
                    end
                catch err
                    if isa(err, ParamBoundsError)
                        draw_loglh = draw_logpost = -Inf
                    #elseif isa(err, SPDError)
                    #    draw_loglh = draw_logpost = -Inf
                    else
                     #   draw_loglh = draw_logpost = -Inf
                        throw(err)
                    end
                end

                if isinf(draw_loglh)
                    draw = vec(rand(m.parameters, 1))
                else
                    success = true
                end
            end
            vector_reshape(draw, draw_loglh, draw_logpost)
        end
        loglh = dropdims(loglh, dims = 1)
        logpost = dropdims(logpost, dims = 1)
    else
        draws = rand(m.parameters, n_parts)
        for i in 1:n_parts
            success = false
            while !success
                try
                    update!(m, draws[:, i])
                    loglh[i] = likelihood(m, data, catch_errors = true, use_chand_recursion = use_chand_recursion, verbose = verbose)
                    logpost[i] = prior(m)
                    if (loglh[i] == -Inf) | (loglh[i]===NaN)
                        logpost[i] = -Inf
                        loglh[i] = -Inf
                    end
                catch err
                    if isa(err, ParamBoundsError)
                        loglh[i] = logpost[i] = -Inf #draws[:, i] = rand(m.parameters, 1)
                        #continue
                    else
                        throw(err)
                    end
                end
                if isinf(loglh[i])
                    draws[:, i] = rand(m.parameters, 1)
                else
                    success = true
                end
            end
        end
    end

    update_draws!(c, draws)
    update_loglh!(c, loglh)
    update_logpost!(c, logpost)
end

# This function is made for transfering the log-likelihood values saved in the
# ParticleCloud from a previous estimation to each particle's respective old_loglh
# field, and for evaluating/saving the likelihood and posterior at the new data, which
# here is just the argument, data.
function initialize_likelihoods!(m::AbstractModel, data::Matrix{Union{Float64, Missing}}, c::ParticleCloud;
                                 parallel::Bool = false, verbose::Symbol = :low)
    # Retire log-likelihood values from the old estimation to the field old_loglh
    map(p -> p.old_loglh = p.loglh, c.particles)

    n_parts = length(c)
    draws = get_vals(c)
    loglh = zeros(n_parts)
    logpost = zeros(n_parts)

    if parallel
        loglh, logpost = @distributed (scalar_reduce) for i in 1:n_parts
            update!(m, draws[:, i])
            draw_loglh = likelihood(m, data, verbose = verbose)
            draw_logpost = prior(m)
            scalar_reshape(draw_loglh, draw_logpost)
        end
    else
        for i in 1:n_parts
            update!(m, draws[:, i])
            loglh[i] = likelihood(m, data, verbose = verbose)
            logpost[i] = prior(m)

            # Will need a way to handle the case when the likelihood with the new data
            # cannot be evaluated (returning -Inf) even if the likelihood was not -Inf
            # prior to incorporating the new data
        end
    end
    update_loglh!(c, loglh)
    update_logpost!(c, logpost)
end

function initialize_cloud_settings!(m::AbstractModel, cloud::ParticleCloud; tempered_update::Bool = false)
    n_parts = length(cloud)

    cloud.tempering_schedule = zeros(1)
    if tempered_update
        cloud.ESS = [cloud.ESS[end]]
    else
        cloud.ESS[1] = n_parts
    end
    cloud.stage_index = 1
    cloud.n_Φ = get_setting(m, :n_Φ)
    cloud.resamples = 0
    cloud.c = get_setting(m, :step_size_smc)
    cloud.accept = get_setting(m, :target_accept)
    cloud.total_sampling_time = 0.

    return nothing
end
