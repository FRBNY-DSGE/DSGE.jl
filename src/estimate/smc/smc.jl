"""
```
smc(m::AbstractDSGEModel, data::Matrix; verbose::Symbol, old_data::Matrix)
smc(m::AbstractDSGEModel, data::DataFrame)
smc(m::AbstractDSGEModel)
```

### Arguments:

- `m`: A model object, which stores parameter values, prior dists, bounds, and various
    other settings that will be referenced
- `data`: A matrix or dataframe containing the time series of the observables used in
    the calculation of the posterior/likelihood
- `old_data`: A matrix containing the time series of observables of previous data
    (with `data` being the new data) for the purposes of a time tempered estimation
    (that is, using the posterior draws from a previous estimation as the initial set
    of draws for an estimation with new data)

### Keyword Arguments:
- `verbose`: Desired frequency of function progress messages printed to standard out.
	- `:none`: No status updates will be reported.
	- `:low`: Status updates for SMC initialization and recursion will be included.
	- `:high`: Status updates for every iteration of SMC is output, which includes
    the mean and standard deviation of each parameter draw after each iteration,
    as well as calculated acceptance rate, ESS, and number of times resampled.
- `log_prob_oldy::Float64 = 0.0`: Log p(\tilde y) which is the log marginal data density
    of the bridge estimation.

### Outputs

- `cloud`: The ParticleCloud object containing all of the information about the
    parameter values from the sample, their respective log-likelihoods, the ESS
    schedule, tempering schedule etc., which is saved in the saveroot.

### Overview

Sequential Monte Carlo can be used in lieu of Random Walk Metropolis Hastings to
    generate parameter samples from high-dimensional parameter spaces using
    sequentially constructed proposal densities to be used in iterative importance
    sampling.

The implementation here is based on Edward Herbst and Frank Schorfheide's 2014 paper
    'Sequential Monte Carlo Sampling for DSGE Models' and the code accompanying their
    book 'Bayesian Estimation of DSGE Models'.

SMC is broken up into three main steps:

- `Correction`: Reweight the particles from stage n-1 by defining incremental weights,
    which gradually "temper in" the likelihood function p(Y|θ)^(ϕ_n - ϕ_n-1) into the
    normalized particle weights.
- `Selection`: Resample the particles if the distribution of particles begins to
    degenerate, according to a tolerance level for the ESS.
- `Mutation`: Propagate particles {θ(i), W(n)} via N(MH) steps of a Metropolis
    Hastings algorithm.
"""
function smc(m::AbstractDSGEModel, data::Matrix{Float64};
             verbose::Symbol = :low,
             old_data::Matrix{Float64} = Matrix{Float64}(undef, size(data, 1), 0),
             old_cloud::ParticleCloud  = ParticleCloud(m, 0),
             run_test::Bool = false,
             filestring_addl::Vector{String} = Vector{String}(),
             continue_intermediate::Bool = false, intermediate_stage_start::Int = 0,
             save_intermediate::Bool = false, intermediate_stage_increment::Int = 10,
             log_prob_oldy::Float64 = 0.0)

    ########################################################################################
    ### Setting Parameters
    ########################################################################################
    # Construct closure of mutation function so as to avoid issues with serialization
    # across workers with different Julia system images
    sendto(workers(), m = m)
    sendto(workers(), data = data)
    function mutation_closure(p::Vector{S}, d_μ::Vector{S}, d_Σ::Matrix{S},
                              blocks_free::Vector{Vector{Int64}},
                              blocks_all::Vector{Vector{Int64}},
                              ϕ_n::S, ϕ_n1::S; c::S = 1.0, α::S = 1.0,
                              old_data::T = Matrix{S}(undef, size(data, 1), 0),
                              use_chand_recursion::Bool = false,
                              verbose::Symbol = :low) where {S<:Float64, T<:AbstractMatrix}
        return mutation(m, data, p, d_μ, d_Σ, blocks_free, blocks_all, ϕ_n, ϕ_n1; c = c, α = α,
                        old_data = old_data, use_chand_recursion = use_chand_recursion,
                        verbose = verbose)
    end
    @everywhere function mutation_closure(p::Vector{S}, d_μ::Vector{S}, d_Σ::Matrix{S},
                                          blocks_free::Vector{Vector{Int64}},
                                          blocks_all::Vector{Vector{Int64}},
                                          ϕ_n::S, ϕ_n1::S; c::S = 1.0, α::S = 1.0,
                                          old_data::T = Matrix{S}(undef, size(data, 1), 0),
                                          use_chand_recursion::Bool = false,
                                          verbose::Symbol = :low) where {S<:Float64, T<:Matrix}
        return mutation(m, data, p, d_μ, d_Σ, blocks_free, blocks_all, ϕ_n, ϕ_n1; c = c, α = α,
                        old_data = old_data, use_chand_recursion = use_chand_recursion,
                        verbose = verbose)
    end

    # General
    parallel = get_setting(m, :use_parallel_workers)
    n_parts  = get_setting(m, :n_particles)
    n_blocks = get_setting(m, :n_smc_blocks)
    n_steps  = get_setting(m, :n_mh_steps_smc)
    n_params = n_parameters(m)

    use_chand_recursion = get_setting(m, :use_chand_recursion)
    if any(isnan.(data)) & use_chand_recursion
        error("Cannot use Chandrasekhar recursions with missing data")
    end

    # Time Tempering
    tempered_update = !isempty(old_data)

    # Check that if there's a tempered update, old and current vintages are different
    if tempered_update
        old_vintage = get_setting(m, :previous_data_vintage)
        @assert old_vintage != data_vintage(m)
    end

    # Step 1 (ϕ schedule) settings
    i = 1                         # Index tracking the stage of the algorithm
    j = 2                         # Index tracking the fixed_schedule entry ϕ_prop
    ϕ_n    = 0.                   # Instantiate ϕ_n and ϕ_prop variables for
    ϕ_prop = 0.                   # reference in respective while loop conditions
    resampled_last_period = false # Ensures proper resetting of ESS_bar after resample
    λ      = get_setting(m, :λ)
    n_Φ    = get_setting(m, :n_Φ)
    tempering_target   = get_setting(m, :adaptive_tempering_target_smc)
    use_fixed_schedule = tempering_target == 0.0

    # Step 2 (Correction) settings
    resampling_method = get_setting(m, :resampler_smc)
    threshold_ratio   = get_setting(m, :resampling_threshold)
    threshold         = threshold_ratio * n_parts

    # Step 3 (Mutation) settings
    c      = get_setting(m, :step_size_smc)
    α      = get_setting(m, :mixture_proportion)
    target = accept = get_setting(m, :target_accept)

    para_symbols    = [θ.key for θ in m.parameters]
    fixed_para_inds = findall([θ.fixed for θ in m.parameters])
    free_para_inds  = findall([!θ.fixed for θ in m.parameters])
    n_free_para     = length(free_para_inds)

    # Step 4 (Initialization of Particle Array Cloud)
    cloud = Cloud(m, n_parts)

    #################################################################################
    ### Initialize Algorithm: Draws from prior
    #################################################################################

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("\n\n SMC " * (m.testing ? "testing " : "") * "starts ....\n\n")
    end

    if tempered_update
        cloud = if isempty(old_cloud)
            loadpath = rawpath(m, "estimate", "smc_cloud.jld2", filestring_addl)
            loadpath = replace(loadpath, r"vint=[0-9]{6}", "vint=" * old_vintage)
            Cloud(load(loadpath, "cloud"))
        else
            Cloud(old_cloud)
        end
        initialize_cloud_settings!(m, cloud; tempered_update = tempered_update)
        initialize_likelihoods!(m, data, cloud, parallel = parallel, verbose = verbose)
    elseif continue_intermediate
        loadpath = rawpath(m, "estimate", "smc_cloud_stage=$(intermediate_stage_start).jld2",
                           filestring_addl)
        cloud    = Cloud(load(loadpath, "cloud"))
    else
        # Instantiating Cloud object, update draws, loglh, & logpost
        initial_draw!(m, data, cloud, parallel = parallel,
                      use_chand_recursion = use_chand_recursion, verbose = verbose)
        initialize_cloud_settings!(m, cloud; tempered_update = tempered_update)
    end

    # Fixed schedule for construction of ϕ_prop
    if use_fixed_schedule
        cloud.tempering_schedule = ((collect(1:n_Φ) .- 1) / (n_Φ-1)) .^ λ
    else
        proposed_fixed_schedule  = ((collect(1:n_Φ) .- 1) / (n_Φ-1)) .^ λ
    end

    # Instantiate incremental and normalized weight matrices for logMDD calculation
    if continue_intermediate
        z_matrix = load(loadpath, "z")
        w_matrix = load(loadpath, "w")
        W_matrix = load(loadpath, "W")
        j        = load(loadpath, "j")
        i        = cloud.stage_index
        c        = cloud.c
        ϕ_prop   = proposed_fixed_schedule[j]
    else
        z_matrix = ones(1)
        w_matrix = zeros(n_parts, 1)
        W_matrix = tempered_update ? get_weights(cloud) : fill(1/n_parts, (n_parts, 1))
    end

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        init_stage_print(cloud, para_symbols; verbose = verbose,
                         use_fixed_schedule = use_fixed_schedule)
    end

    #################################################################################
    ### Recursion
    #################################################################################
    (VERBOSITY[verbose] >= VERBOSITY[:low]) && println("\n\n SMC recursion starts \n\n")

    while ϕ_n < 1.

        start_time = time_ns()
        cloud.stage_index = i += 1

        #############################################################################
        ### Step 0: Setting ϕ_n (either adaptively or by the fixed schedule)
        #############################################################################
        ϕ_n1 = cloud.tempering_schedule[i-1]

        if use_fixed_schedule
            ϕ_n = cloud.tempering_schedule[i]
        else
            ϕ_n, resampled_last_period, j, ϕ_prop = solve_adaptive_ϕ(cloud,
                                                       proposed_fixed_schedule,
                                                       i, j, ϕ_prop, ϕ_n1,
                                                       tempering_target,
                                                       resampled_last_period)
        end

        #############################################################################
        ### Step 1: Correction
        #############################################################################
        # Calculate incremental weights (if no old data, get_old_loglh(cloud) = 0)
        if tempered_update_prior_weight == 0.0
            incremental_weights = exp.((ϕ_n1 - ϕ_n) * get_old_loglh(cloud) +
                                       (ϕ_n - ϕ_n1) * get_loglh(cloud))
        elseif tempered_update_prior_weight == 1.0
            incremental_weights = exp.((ϕ_n - ϕ_n1) * get_loglh(cloud))
        else
            incremental_weights = exp.((ϕ_n1 - ϕ_n) * log.((exp.(get_old_loglh(cloud) .- log_prob_oldy .+
                                                       log(1-tempered_update_prior_weight)) .+ tempered_update_prior_weight)) +
                                       (ϕ_n - ϕ_n1) * get_loglh(cloud))
        end

        # Update weights
        update_weights!(cloud, incremental_weights)
        mult_weights = get_weights(cloud)

        # Normalize weights
        normalize_weights!(cloud)
        normalized_weights = get_weights(cloud)

        push!(z_matrix, sum(mult_weights))
        w_matrix = hcat(w_matrix, incremental_weights)
        W_matrix = hcat(W_matrix, normalized_weights)

        ##############################################################################
        ### Step 2: Selection
        ##############################################################################

        # Calculate the degeneracy/effective sample size metric
        push!(cloud.ESS, 1 / sum(normalized_weights .^ 2))

        # If this assertion does not hold then there are likely too few particles
        @assert !isnan(cloud.ESS[i]) "no particles have non-zero weight"

        # Resample if degeneracy/ESS metric falls below the accepted threshold
        if (cloud.ESS[i] < threshold)
            # Resample according to particle weights, uniformly reset weights to 1/n_parts
            new_inds = resample(normalized_weights; method = resampling_method)
            cloud.particles = [deepcopy(cloud.particles[i,j]) for i in new_inds,
                               j=1:size(cloud.particles, 2)]
            reset_weights!(cloud)
            cloud.resamples += 1
            resampled_last_period = true
            W_matrix[:, i] = fill(1/n_parts, (n_parts, 1))
        end

        ##############################################################################
        ### Step 3: Mutation
        ##############################################################################

        # Calculate adaptive c-step for use as scaling coefficient in mutation MH step
        c = c * (0.95 + 0.10 * exp(16.0 * (cloud.accept - target)) /
                 (1.0 + exp(16.0 * (cloud.accept - target))))
        cloud.c = c

        θ_bar = weighted_mean(cloud)
        R     = weighted_cov(cloud)

        # Ensures marix is positive semi-definite symmetric
        # (not off due to numerical error) and values haven't changed
        R_fr = (R[free_para_inds, free_para_inds] + R[free_para_inds, free_para_inds]') / 2.

        # MvNormal centered at ̄θ with var-cov ̄Σ, subsetting out the fixed parameters
        θ_bar_fr = θ_bar[free_para_inds]

        # Generate random parameter blocks
        blocks_free = generate_free_blocks(n_free_para, n_blocks)
        blocks_all  = generate_all_blocks(blocks_free, free_para_inds)

        new_particles = if parallel
            @distributed (hcat) for k in 1:n_parts
                mutation_closure(cloud.particles[k, :], θ_bar_fr, R_fr, blocks_free,
                                 blocks_all, ϕ_n, ϕ_n1; c = c, α = α, old_data = old_data,
                                 use_chand_recursion = use_chand_recursion,
                                 verbose = verbose)
            end
        else
            hcat([mutation_closure(cloud.particles[k, :], θ_bar_fr, R_fr,
                                   blocks_free, blocks_all, ϕ_n, ϕ_n1; c = c,
                                   α = α, old_data = old_data,
                                   use_chand_recursion = use_chand_recursion,
                                   verbose = verbose) for k=1:n_parts]...)
        end
        update_cloud!(cloud, new_particles)
        update_acceptance_rate!(cloud)

        ##############################################################################
        ### Timekeeping and Output Generation
        ##############################################################################
        cloud.total_sampling_time += Float64((time_ns() - start_time) * 1e-9)

        if VERBOSITY[verbose] >= VERBOSITY[:low]
            end_stage_print(cloud, para_symbols; verbose = verbose,
                            use_fixed_schedule = use_fixed_schedule)
        end

        if run_test && (i == 3)
            break
        end
        if mod(cloud.stage_index, intermediate_stage_increment) == 0 && save_intermediate
            jldopen(rawpath(m, "estimate", "smc_cloud_stage=$(cloud.stage_index).jld2"),
                    true, true, true, IOStream) do file
                write(file, "cloud", ParticleCloud(cloud, para_symbols))
                write(file, "w", w_matrix)
                write(file, "W", W_matrix)
                write(file, "z", z_matrix)
                write(file, "j", j)
            end
        end
    end

    ##################################################################################
    ### Saving data
    ##################################################################################

    if !m.testing
        simfile = h5open(rawpath(m, "estimate", "smcsave.h5", filestring_addl), "w")
        particle_store = d_create(simfile, "smcparams", datatype(Float64),
                                  dataspace(n_parts, n_params))
        for i in 1:n_parts; particle_store[i,:] = cloud.particles[i, 1:n_params] end
        close(simfile)

        jldopen(rawpath(m, "estimate", "smc_cloud.jld2", filestring_addl),
                true, true, true, IOStream) do file
            write(file, "cloud", ParticleCloud(cloud, para_symbols))
            write(file, "w", w_matrix)
            write(file, "W", W_matrix)
            write(file, "z", z_matrix)
        end
    end
end

function smc(m::AbstractDSGEModel, data::DataFrame; verbose::Symbol = :low,
             save_intermediate::Bool = false, intermediate_stage_increment::Int = 10,
             filestring_addl::Vector{String} = Vector{String}(undef, 0), log_prob_oldy::Float64 = 0.0)
    data_mat = df_to_matrix(m, data)
    return smc(m, data_mat, verbose = verbose, save_intermediate = save_intermediate,
               filestring_addl = filestring_addl, log_prob_oldy = log_prob_oldy)
end

function smc(m::AbstractDSGEModel; verbose::Symbol = :low,
             save_intermediate::Bool = false, intermediate_stage_increment::Int = 10,
             filestring_addl::Vector{String} = Vector{String}(undef, 0), log_prob_oldy::Float64 = 0.0)
    data = load_data(m)
    data_mat = df_to_matrix(m, data)
    return smc(m, data_mat, verbose=verbose, save_intermediate = save_intermediate,
               filestring_addl = filestring_addl, log_prob_oldy = log_prob_oldy)
end
