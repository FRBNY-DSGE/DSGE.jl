"""
```
initial_draw(m::AbstractModel, data::Matrix{Float64}, c::ParticleCloud)
```

Draw from a general starting distribution (set by default to be from the prior) to initialize the SMC algorithm.
Returns a tuple (logpost, loglh) and modifies the particle objects in the particle cloud in place.

"""
function initial_draw(m::AbstractModel, data::Matrix{Float64}, c::ParticleCloud;
                      parallel::Bool = false)
    n_parts = length(c)
    loglh = zeros(n_parts)
    logpost = zeros(n_parts)
    if parallel
        out = @parallel (hcat) for i in 1:n_parts
            draw = vec(rand(m.parameters, 1))
            draw_loglh = 0.
            draw_logpost = 0.
            success = false
            while !success
                try
                    update!(m, draw)
                    draw_loglh   = likelihood(m, data)
                    draw_logpost = prior(m)
                catch
                    draw = vec(rand(m.parameters, 1))
                    continue # Keep drawing until you get valid draws
                end
                success = true
            end
            (draw, draw_loglh, draw_logpost)
        end
        draws = zeros(n_parameters(m), n_parts)
        for i in 1:n_parts
            draws[:, i] = out[i][1]
            loglh[i]    = out[i][2]
            logpost[i]  = out[i][3]
        end
    else
        draws = rand(m.parameters, n_parts)
        for i in 1:n_parts
            success = false
            while !success
                try
                    update!(m, draws[:, i])
                    loglh[i] = likelihood(m, data)
                    logpost[i] = prior(m)
                catch
                    draws[:, i] = rand(m.parameters, 1)
                    continue
                end
                success = true
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
function initialize_likelihoods(m::AbstractModel, data::Matrix{Float64}, c::ParticleCloud;
                                parallel::Bool = false)
    # Retire log-likelihood values from the old estimation to the field old_loglh
    map(p -> p.old_loglh = p.loglh, c.particles)

    n_parts = length(c)
    draws = get_vals(c)
    loglh = zeros(n_parts)
    logpost = zeros(n_parts)

    if parallel
        out = @sync @parallel (hcat) for i in 1:n_parts
            update!(m, draws[:, i])
            draw_loglh = likelihood(m, data)
            draw_logpost = prior(m)
            (draw_loglh, draw_logpost)
        end
        for i in 1:n_parts
            loglh[i]    = out[i][1]
            logpost[i]  = out[i][2]
        end
    else
        for i in 1:n_parts
            update!(m, draws[:, i])
            loglh[i] = likelihood(m, data)
            logpost[i] = prior(m)

            # Will need a way to handle the case when the likelihood with the new data
            # cannot be evaluated (returning -Inf) even if the likelihood was not -Inf
            # prior to incorporating the new data
        end
    end
    update_loglh!(c, loglh)
    update_logpost!(c, logpost)
end

function solve_adaptive_ϕ(cloud::ParticleCloud, proposed_fixed_schedule::Vector{Float64},
                          i::Int64, j::Int64, ϕ_prop::Float64, ϕ_n1::Float64,
                          tempering_target::Float64, n_Φ::Int64, endo_type::Symbol,
                          resampled_last_period::Bool)
    if resampled_last_period
        # The ESS_bar should be reset to target an evenly weighted particle population
        ESS_bar = tempering_target*length(cloud)
        resampled_last_period = false
    else
        ESS_bar = tempering_target*cloud.ESS[i-1]
    end

    # Setting up the optimal ϕ solving function for endogenizing the tempering schedule
    optimal_ϕ_function(ϕ)  = compute_ESS(get_loglh(cloud), get_weights(cloud), ϕ, ϕ_n1,
                                         old_loglh = get_old_loglh(cloud)) - ESS_bar

    # Find a ϕ_prop such that the optimal ϕ_n will lie somewhere between ϕ_n1 and ϕ_prop
    # Do so by iterating through a proposed_fixed_schedule and finding the first
    # ϕ_prop such that the ESS would fall by more than the targeted amount ESS_bar
    if endo_type == :adaptive
        while optimal_ϕ_function(ϕ_prop) >= 0 && j <= n_Φ
            ϕ_prop = proposed_fixed_schedule[j]
            j += 1
        end
    elseif endo_type == :adaptive_min
        ϕ_prop = proposed_fixed_schedule[j]
    elseif endo_type == :adaptive_min_j_eq_n
        if j <= n_Φ
            ϕ_prop = proposed_fixed_schedule[j]
        else
            ϕ_prop = 1.
        end
    else
        throw("not an endo type")
    end

    # Note: optimal_ϕ_function(ϕ_n1) > 0 because ESS_{t-1} is always positive
    # When ϕ_prop != 1. then there are still ϕ increments strictly below 1 that
    # give the optimal ϕ step, ϕ_n.
    # When ϕ_prop == 1. but optimal_ϕ_function(ϕ_prop) < 0 then there still exists
    # an optimal ϕ step, ϕ_n, that does not equal 1.
    # Thus the interval [optimal_ϕ_function(ϕ_n1), optimal_ϕ_function(ϕ_prop)] always
    # contains a 0 by construction.

    # Modification makes it such that ϕ_n is the minimum of ϕ_prop (the fixed schedule)
    # at a given stage or the root-solved ϕ such that the ESS drops by the target amount.
    # Thus the ϕ_schedule should be strictly bounded above by the fixed schedule
    # i.e. the adaptive ϕ schedule should not outpace the fixed schedule at the end
    # (when the fixed schedule tends to drop by less than 5% per iteration)

    if ϕ_prop != 1. || optimal_ϕ_function(ϕ_prop) < 0
        if endo_type == :adaptive
            ϕ_n = fzero(optimal_ϕ_function, [ϕ_n1, ϕ_prop], xtol = 0.)
        elseif endo_type == :adaptive_min
            if optimal_ϕ_function(ϕ_prop) < 0
                ϕ_n = fzero(optimal_ϕ_function, [ϕ_n1, ϕ_prop], xtol = 0.)
            else
                ϕ_n = ϕ_prop
                j += 1
            end
        elseif endo_type == :adaptive_min_j_eq_n
            if optimal_ϕ_function(ϕ_prop) < 0
                ϕ_n = fzero(optimal_ϕ_function, [ϕ_n1, ϕ_prop], xtol = 0.)
            else
                ϕ_n = ϕ_prop
            end
            j += 1
        else
            throw("not an endo type")
        end
        push!(cloud.tempering_schedule, ϕ_n)
    else
        ϕ_n = 1.
        push!(cloud.tempering_schedule, ϕ_n)
    end

    return ϕ_n, resampled_last_period, j, ϕ_prop
end

"""
```
mvnormal_mixture_draw{T<:AbstractFloat}(θ_old, σ; cc, α, θ_prop)
```

Create a `DegenerateMvNormal` distribution object, `d`, from a parameter vector, `p`, and a
standard deviation matrix (obtained from SVD), `σ`.

Generate a draw from the mixture distribution of:
1. A `DegenerateMvNormal` centered at θ_old with the standard deviation matrix `σ`, scaled by `cc^2` and with mixture proportion `α`.
2. A `DegenerateMvNormal` centered at the same mean, but with a standard deviation matrix of the diagonal entries of `σ` scaled by `cc^2` with mixture proportion `(1 - α)/2`.
3. A `DegenerateMvNormal`  with the same standard deviation matrix `σ` but centered at the new proposed mean, `θ_prop`, scaled by `cc^2`, and with mixture proportion `(1 - α)/2`.

If no `θ_prop` is given, but an `α` is specified, then the mixture will consist of `α` of
the standard distribution and `(1 - α)` of the diagonalized distribution.

### Arguments
- `θ_old::Vector{T}`: The mean of the desired distribution
- `σ::Matrix{T}`: The standard deviation matrix of the desired distribution

### Keyword Arguments
- `cc::T`: The standard deviation matrix scaling factor
- `α::T`: The mixing proportion
- `θ_prop::Vector{T}`: The proposed parameter vector to be used as part of the mixture distribution, set by default to be the weighted mean of the particles, prior to mutation.

### Outputs
- `θ_new::Vector{T}`: The draw from the mixture distribution to be used as the MH proposed step
- `new mixture_density::T`: The mixture density conditional on θ_old evaluated at `θ_new` to be used in calculating the MH move probability
- `old mixture_density::T`: The mixture density conditional on θ_new evaluated at `θ_old` to be used in calculating the MH move probability

"""
function mvnormal_mixture_draw{T<:AbstractFloat}(θ_old::Vector{T}, d_prop::Distribution;
                                                 cc::T = 1.0, α::T = 1.)
    @assert 0 <= α <= 1

    # Create mixture distribution conditional on the previous parameter value, θ_old
    d_old = MvNormal(θ_old, cc*d_prop.Σ)
    d_diag_old = MvNormal(θ_old, diagm(diag(cc*d_prop.Σ)))
    d_mix_old = MixtureModel(MvNormal[d_old, d_diag_old, d_prop], [α, (1 - α)/2, (1 - α)/2])

    θ_new = rand(d_mix_old)

    # Create mixture distribution conditional on the new parameter value, θ_new
    d_new = MvNormal(θ_new, cc*d_prop.Σ)
    d_diag_new = MvNormal(θ_new, diagm(diag(cc*d_prop.Σ)))
    d_mix_new = MixtureModel(MvNormal[d_new, d_diag_new, d_prop], [α, (1 - α)/2, (1 - α)/2])

    # To clarify, this is not just the density of θ_new/θ_old using a given mixture
    # density function, but rather, the density of θ_new | θ_old and the density of
    # θ_old | θ_new taken with respect to their respective mixture densities
    new_mixture_density = logpdf(d_mix_old, θ_new)
    old_mixture_density = logpdf(d_mix_new, θ_old)

    return θ_new, new_mixture_density, old_mixture_density
end

function compute_ESS{T<:AbstractFloat}(loglh::Vector{T}, current_weights::Vector{T},
                                       ϕ_n::T, ϕ_n1::T;
                                       old_loglh::Vector{T} = zeros(length(loglh)))
    incremental_weights = exp.((ϕ_n1 - ϕ_n)*old_loglh + (ϕ_n - ϕ_n1)*loglh)
    new_weights = current_weights.*incremental_weights
    normalized_weights = new_weights/sum(new_weights)
    ESS = 1/sum(normalized_weights.^2)
    return ESS
end

# For mutation
# Generate a Vector of Vector{Int64} of length n_blocks, where each
# element contains a subset of the randomly permuted set of indices 1:n_para
# For the purpose of indexing into R_fr
function generate_free_blocks(n_para::Int64, n_blocks::Int64)
    rand_inds = shuffle(1:n_para)

    subset_length = cld(n_para, n_blocks) # ceiling division
    last_block_length = n_para - subset_length*(n_blocks - 1)

    blocks_free = Vector{Vector{Int64}}(n_blocks)
    for i in 1:n_blocks
        if i < n_blocks
            blocks_free[i] = rand_inds[((i-1)*subset_length + 1):(i*subset_length)]
        else
            # To account for the fact that the last block may be smaller than the others
            blocks_free[i] = rand_inds[end-last_block_length+1:end]
        end
    end
    return blocks_free
end

# For mutation
# Generate a Vector of Vector{Int64} of length n_blocks, where each
# element contains a subset corresponding to the subset of blocks_free of the same
# index but with indices that map to free_para_inds as opposed to 1:n_para
# For the purpose of "re-creating" the proposed parameter vector that contains both free
# and fixed parameters from the mh step generated from only the free parameters
function generate_all_blocks(blocks_free::Vector{Vector{Int64}}, free_para_inds::Vector{Int64})
    n_free_para = length(free_para_inds)
    # Need to know the mapping from an ordered list of 1:n_free_para
    # to the index in the actual parameter vector
    ind_mappings = Dict{Int64, Int64}()
    for (k, v) in zip(1:n_free_para, free_para_inds)
        ind_mappings[k] = v
    end

    # Want: Input: blocks, a vector of vectors of indices of randomized blocks of an ordered list of 1:n_free_para
    # Output: rev_blocks, a vector of vector of indices of the rand blocks indices that correspond to the actual
    # indices in a parameter vector
    function block_map(blocks::Vector{Vector{Int64}}, ind_mappings::Dict{Int64, Int64})
        blocks_all = similar(blocks)
        for (i, block) in enumerate(blocks)
            blocks_all[i] = similar(block)
            for (j, b) in enumerate(block)
                blocks_all[i][j] = ind_mappings[b]
            end
        end
        return blocks_all
    end

    blocks_all = block_map(blocks_free, ind_mappings)

    return blocks_all
end

function init_stage_print(cloud::ParticleCloud;
                          verbose::Symbol=:low, use_fixed_schedule::Bool = true)
    if use_fixed_schedule
        println("--------------------------")
            println("Iteration = $(cloud.stage_index) / $(cloud.n_Φ)")
    else
        println("--------------------------")
            println("Iteration = $(cloud.stage_index)")
    end
	println("--------------------------")
        println("phi = $(cloud.tempering_schedule[cloud.stage_index])")
	println("--------------------------")
        println("c = $(cloud.c)")
        println("ESS = $(cloud.ESS[cloud.stage_index])   ($(cloud.resamples) total resamples.)")
	println("--------------------------")
    if VERBOSITY[verbose] >= VERBOSITY[:high]
        μ = weighted_mean(cloud)
        σ = weighted_std(cloud)
        for n=1:length(cloud.particles[1])
            println("$(cloud.particles[1].keys[n]) = $(round(μ[n], 5)), $(round(σ[n], 5))")
	    end
    end
end

function end_stage_print(cloud::ParticleCloud;
                         verbose::Symbol=:low, use_fixed_schedule::Bool = true)
    total_sampling_time_minutes = cloud.total_sampling_time/60
    if use_fixed_schedule
        expected_time_remaining_sec = (cloud.total_sampling_time/cloud.stage_index)*(cloud.n_Φ - cloud.stage_index)
        expected_time_remaining_minutes = expected_time_remaining_sec/60
    end

    println("--------------------------")
    if use_fixed_schedule
        println("Iteration = $(cloud.stage_index) / $(cloud.n_Φ)")
        println("time elapsed: $(round(total_sampling_time_minutes, 4)) minutes")
        println("estimated time remaining: $(round(expected_time_remaining_minutes, 4)) minutes")
    else
        println("Iteration = $(cloud.stage_index)")
        println("time elapsed: $(round(total_sampling_time_minutes, 4)) minutes")
    end
    println("--------------------------")
        println("phi = $(cloud.tempering_schedule[cloud.stage_index])")
    println("--------------------------")
        println("c = $(cloud.c)")
        println("accept = $(cloud.accept)")
        println("ESS = $(cloud.ESS[cloud.stage_index])   ($(cloud.resamples) total resamples.)")
    println("--------------------------")
    if VERBOSITY[verbose] >= VERBOSITY[:high]
        μ = weighted_mean(cloud)
        σ = weighted_std(cloud)
        for n=1:length(cloud.particles[1])
            println("$(cloud.particles[1].keys[n]) = $(round(μ[n], 5)), $(round(σ[n], 5))")
        end
    end
end
