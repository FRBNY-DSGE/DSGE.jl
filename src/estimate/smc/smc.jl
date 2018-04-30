"""
```
smc(m::AbstractModel, data::Matrix; verbose::Symbol, old_data::Matrix)
smc(m::AbstractModel, data::DataFrame)
smc(m::AbstractModel)
```

### Arguments:

- `m`: A model object, from which its parameters values, prior dists, and bounds will be referenced
- `data`: A matrix or data frame containing the time series of the observables to be used in the calculation of the posterior/likelihood
- `old_data`: A matrix containing the time series of observables of previous data (with `data` being the new data) for the purposes of time tempered estimation

### Keyword Arguments:
- `verbose`: The desired frequency of function progress messages printed to standard out.
	- `:none`: No status updates will be reported.
	- `:low`: Status updates for SMC initialization and recursion will be included.
	- `:high`: Status updates for every iteration of SMC is output, which includes the mu and sigma of each individual parameter draw after each iteration, as well as the calculated acceptance rate, the ESS, and the number of times resampled.

### Outputs

- `cloud`: The ParticleCloud object that contains all of the information about the parameter values from the sample, their respective log-likelihoods, the ESS schedule, tempering schedule etc. is saved in the saveroot.

### Overview

Sequential Monte Carlo can be used in lieu of Random Walk Metropolis Hastings to generate parameter samples from high-dimensional parameter spaces using sequentially constructed proposal densities.

The implementation here is based upon Edward Herbst and Frank Schorfheide's 2014 paper 'Sequential Monte Carlo Sampling for DSGE Models' and the code accompanying their book 'Bayesian Estimation of DSGE Models'.

SMC is broken up into three main steps:

- `Correction`: Reweight the particles from stage n-1 by defining "incremental weights", incweight, which gradually incorporate the likelihood function p(Y|θ(i, n-1)) into the particle weights.
- `Selection`: Resample the particles if the distribution of particles begins to degenerate, according to a tolerance level ESS < N/2. The current resampling technique employed is systematic resampling.

- `Mutation`: Propagate particles {θ(i), W(n)} via N(MH) steps of a Metropolis Hastings algorithm.
"""
function smc(m::AbstractModel, data::Matrix{Float64};
             verbose::Symbol = :low, old_data::Matrix{Float64} = Matrix{Float64}(size(data, 1), 0))
    ########################################################################################
    ### Setting Parameters
    ########################################################################################

    # Temporary
    endo_type = get_setting(m, :endogenous_type)

    # General
    parallel = get_setting(m, :use_parallel_workers)
    n_parts = get_setting(m, :n_particles)
    n_params = n_parameters(m)
    n_blocks = get_setting(m, :n_smc_blocks)

    # Time Tempering
    tempered_update = !isempty(old_data)
    # Quick check that if there is a tempered update that the old vintage and current vintage are different
    if tempered_update
        old_vintage = get_setting(m, :previous_data_vintage)
        @assert old_vintage != data_vintage(m)
    end

    # Step 0 (ϕ schedule) settings
    use_fixed_schedule = get_setting(m, :use_fixed_schedule)
    λ = get_setting(m, :λ)
    n_Φ = get_setting(m, :n_Φ)
    tempering_target = get_setting(m, :tempering_target)

    # Step 2 (Correction) settings
    resampling_method = get_setting(m, :resampler_smc)
    threshold_ratio = get_setting(m, :resampling_threshold)
    threshold = threshold_ratio * n_parts
    use_CESS = get_setting(m, :use_CESS)

    # Step 3 (Mutation) settings
    c = get_setting(m, :step_size_smc)
    accept = get_setting(m, :init_accept)
    target = get_setting(m, :target_accept)
    α = get_setting(m, :mixture_proportion)
    fixed_para_inds = find([θ.fixed for θ in m.parameters])
    free_para_inds = find([!θ.fixed for θ in m.parameters])

    ########################################################################################
    ### Initialize Algorithm: Draws from prior
    ########################################################################################

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        if m.testing
            println("\n\n SMC testing starts ....  \n\n  ")
        else
            println("\n\n SMC starts ....  \n\n  ")
        end
    end

    # fix seed if testing
    if m.testing
        srand(42)
    end

    if tempered_update
        # Load the previous ParticleCloud as the starting point for time tempering
        loadpath = rawpath(m, "estimate", "smc_cloud.jld")
        loadpath = replace(loadpath, r"vint=[0-9]{6}", "vint="*old_vintage)
        loadpath = replace(loadpath, "temp=true", "temp=false")

        cloud = load(loadpath, "cloud")
        reset_cloud_settings!(cloud)
        initialize_likelihoods(m, data, cloud, parallel = parallel)
    else
        # Instantiating ParticleCloud object
        cloud = ParticleCloud(m, n_parts)

        # Modifies the cloud object in place to update draws, loglh, & logpost
        initial_draw(m, data, cloud, parallel = parallel)
    end

    # Fixed schedule for construction of ϕ_prop
    if use_fixed_schedule
        cloud.tempering_schedule = ((collect(1:n_Φ)-1)/(n_Φ-1)).^λ
    else
        proposed_fixed_schedule = ((collect(1:n_Φ)-1)/(n_Φ-1)).^λ
    end

    # Saving initial parameter values in the ParticleCloud instance
    cloud.n_Φ = n_Φ
    cloud.c = c
    cloud.accept = accept
    cloud.stage_index = i = 1
    cloud.ESS[1] = n_parts          # To make adaptive ϕ schedule calculate ESS_bar properly

    j = 1                           # The index tracking to the fixed_schedule entry that ϕ_prop is set as
    resampled_last_period = false   # To ensure proper resetting of ESS_bar right after resample
    ϕ_n = 0.                        # Instantiating ϕ_n and ϕ_prop variables to be referenced in their
    ϕ_prop = 0.                     # respective while loop conditions
    w_matrix = zeros(n_parts, 1)    # Incremental and normalized weight matrices (n_parts x n_Φ) to store
    W_matrix = ones(n_parts, 1)     # for the calculation of the MDD

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        init_stage_print(cloud; verbose = verbose, use_fixed_schedule = use_fixed_schedule)
    end

    ########################################################################################
    ### Recursion
    ########################################################################################

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("\n\n SMC recursion starts \n\n")
    end

    while ϕ_n < 1.

    tic()
    cloud.stage_index = i += 1

    ########################################################################################
    ### Step 0: Setting ϕ_n (either adaptively or by the fixed schedule)
    ########################################################################################
    ϕ_n1 = cloud.tempering_schedule[i-1]

    if use_fixed_schedule
        ϕ_n = cloud.tempering_schedule[i]
    else
        if resampled_last_period
            # The ESS_bar should be reset to target an evenly weighted particle population
            ESS_bar = tempering_target*n_parts
            resampled_last_period = false
        else
            ESS_bar = tempering_target*cloud.ESS[i-1]
        end

        # Setting up the optimal ϕ solving function for endogenizing the tempering schedule
        optimal_ϕ_function(ϕ)  = compute_ESS(get_loglh(cloud), get_weights(cloud), ϕ, ϕ_n1,
                                             use_CESS = use_CESS, old_loglh = get_old_loglh(cloud)) - ESS_bar

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
    end

    ########################################################################################
    ### Step 1: Correction
    ########################################################################################

    # Calculate incremental weights
    incremental_weights = exp.((ϕ_n1 - ϕ_n)*get_old_loglh(cloud) + (ϕ_n - ϕ_n1)*get_loglh(cloud))

    # Need previous weights for CESS calculation
    if use_CESS
        previous_weights = get_weights(cloud)
    end

    # Update weights
    update_weights!(cloud, incremental_weights)

    # Normalize weights
    normalize_weights!(cloud)

    normalized_weights = get_weights(cloud)

    w_matrix = hcat(w_matrix, incremental_weights)
    W_matrix = hcat(W_matrix, normalized_weights)

    ########################################################################################
    ### Step 2: Selection
    ########################################################################################

    # Calculate the degeneracy/effective sample size metric
    if use_CESS
        push!(cloud.ESS, n_parts*sum(previous_weights .* incremental_weights)^2/sum(previous_weights .*
                                                                                    incremental_weights.^2))
    else
        push!(cloud.ESS, 1/sum(normalized_weights.^2))
    end

    # If this assertion does not hold then there are likely too few particles
    @assert !isnan(cloud.ESS[i]) "no particles have non-zero weight"

    # Resample if the degeneracy/effective sample size metric falls below the accepted threshold
    if (cloud.ESS[i] < threshold)
        new_inds = resample(normalized_weights; method = resampling_method,
                            parallel = parallel, testing = m.testing)

        # update parameters/logpost/loglh with resampled values
        # reset the weights to 1/n_parts
        cloud.particles = [deepcopy(cloud.particles[i]) for i in new_inds]
        reset_weights!(cloud)
        cloud.resamples += 1
        resampled_last_period = true
    end

    ########################################################################################
    ### Step 3: Mutation
    ########################################################################################

    # Calculate the adaptive c-step to be used as a scaling coefficient in the mutation MH stea
    c = c*(0.95 + 0.10*exp(16*(cloud.accept - target))/(1 + exp(16*(cloud.accept - target))))
    cloud.c = c

    θ_bar = weighted_mean(cloud)
    R = weighted_cov(cloud)
    R_fr = (R[free_para_inds, free_para_inds] + R[free_para_inds, free_para_inds]')/2

    # MvNormal centered at ̄θ with var-cov ̄Σ, subsetting out the fixed parameters
    d = MvNormal(θ_bar[free_para_inds], R_fr)

    ### BLOCKING ###
    n_para = n_parameters(m)
    n_free_para = length(free_para_inds)
    n_blocks = get_setting(m, :n_smc_blocks)

    # New way of generating blocks
    blocks_free = generate_free_blocks(n_free_para, n_blocks)
    blocks_all  = generate_all_blocks(blocks_free, free_para_inds)

    if parallel
        new_particles = @parallel (vcat) for j in 1:n_parts
            mutation(m, data, cloud.particles[j], d, blocks_free, blocks_all, ϕ_n, ϕ_n1;
                     c = c, α = α, old_data = old_data)
        end
    else
        new_particles = [mutation(m, data, cloud.particles[j], d, blocks_free, blocks_all, ϕ_n, ϕ_n1;
                                  c = c, α = α, old_data = old_data) for j = 1:n_parts]
    end

    cloud.particles = new_particles
    update_acceptance_rate!(cloud) # update average acceptance rate

    ########################################################################################
    ### Timekeeping and Output Generation
    ########################################################################################

    cloud.total_sampling_time += toq()

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        end_stage_print(cloud; verbose = verbose, use_fixed_schedule = use_fixed_schedule)
    end

    end

    ########################################################################################
    ### Saving data
    ########################################################################################

    if !m.testing
        simfile = h5open(rawpath(m, "estimate", "smcsave.h5"),"w")
        particle_store = d_create(simfile, "smcparams", datatype(Float32),
                                  dataspace(n_parts, n_params))
        for i in 1:length(cloud)
            particle_store[i,:] = cloud.particles[i].value
        end
        close(simfile)
        jldopen(rawpath(m, "estimate", "smc_cloud.jld"), "w") do file
            write(file, "cloud", cloud)
            write(file, "w", w_matrix)
            write(file, "W", W_matrix)
        end
    end
end

function smc(m::AbstractModel, data::DataFrame; verbose::Symbol=:low)
    data_mat = df_to_matrix(m, data)
    return smc(m, data_mat, verbose=verbose)
end

function smc(m::AbstractModel; verbose::Symbol=:low)
    data = load_data(m)
    data_mat = df_to_matrix(m, data)
    return smc(m, data_mat, verbose=verbose)
end
