"""
```
smc(m::AbstractModel, data::Matrix)
smc(m::AbstractModel, data::DataFrame)
smc(m::AbstractModel)
smc(m::AbstractModel, cloud::ParticleCloud, old_data::Matrix, new_data::Matrix)
```

### Arguments:

- `m`: A model object, from which its parameters values, prior dists, and bounds will be referenced
- `data`: A matrix or data frame containing time series of the observables to be used in the calculation of the posterior/likelihood

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
function smc(m::AbstractModel, data::Matrix; verbose::Symbol = :low, tempered_update::Bool = false)
    ########################################################################################
    ### Setting Parameters
    ########################################################################################

    parallel = get_setting(m, :use_parallel_workers)
    use_fixed_schedule = get_setting(m, :use_fixed_schedule)

    n_parts = get_setting(m, :n_particles)
    n_params = n_parameters(m)
    n_blocks = get_setting(m, :n_smc_blocks)

    c = get_setting(m, :step_size_smc)
    accept = get_setting(m, :init_accept)
    target = get_setting(m, :target_accept)
    λ = get_setting(m, :λ)
    n_Φ = get_setting(m, :n_Φ)
    tempering_target = get_setting(m, :tempering_target)
    threshold_ratio = get_setting(m, :resampling_threshold)
    threshold = threshold_ratio * n_parts

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

    # if !tempered_update
        # Instantiating ParticleCloud object
        cloud = ParticleCloud(m, n_parts)

        # Particle draws from the parameter's marginal priors
        # Modifies the cloud object in place to update draws, loglh, & logpost
        initial_draw(m, data, cloud)
    # else
        # updated_vintage = get_setting(m, :updated_data_vintage)
        # cloud = load(rawpath(m, "estimate", "smc_cloud_updvint=$updated_vintage.jld"), "old_cloud")

        # unchanged_loglh = load(rawpath(m, "estimate", "smc_cloud_updvint=$updated_vintage.jld"), "unchanged_loglh")
        # revised_loglh = load(rawpath(m, "estimate", "smc_cloud_updvint=$updated_vintage.jld"), "revised_loglh")
        # new_loglh = load(rawpath(m, "estimate", "smc_cloud_updvint=$updated_vintage.jld"), "new_loglh")

        # update_loglh!(cloud, unchanged_loglh)
    # end

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

    total_sampling_time = 0.

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
        optimal_ϕ_function(ϕ)  = compute_ESS(get_loglh(cloud), get_weights(cloud), ϕ, ϕ_n1) - ESS_bar

        # Find a ϕ_prop such that the optimal ϕ_n will lie somewhere between ϕ_n1 and ϕ_prop
        # Do so by iterating through a proposed_fixed_schedule and finding the first
        # ϕ_prop such that the ESS would fall by more than the targeted amount ESS_bar
        while optimal_ϕ_function(ϕ_prop) >= 0 && j <= n_Φ
            ϕ_prop = proposed_fixed_schedule[j]
            j += 1
        end

        # if ϕ_prop == 1.
            # jldopen("debug_file.jld", "w") do file
                # write(file, "cloud", cloud)
                # write(file, "ϕ_prop", ϕ_prop)
                # write(file, "ϕ_n1", ϕ_n1)
                # write(file, "ESS_bar", ESS_bar)
                # write(file, "proposed_fixed_schedule", proposed_fixed_schedule)
            # end
        # end

        # Note: optimal_ϕ_function(ϕ_n1) > 0 because ESS_{t-1} is always positive
        # When ϕ_prop != 1. then there are still ϕ increments strictly below 1 that
        # give the optimal ϕ step, ϕ_n.
        # When ϕ_prop == 1. but optimal_ϕ_function(ϕ_prop) < 0 then there still exists
        # an optimal ϕ step, ϕ_n, that does not equal 1.
        # Thus the interval [optimal_ϕ_function(ϕ_n1), optimal_ϕ_function(ϕ_prop)] always
        # contains a 0 by construction.
        if ϕ_prop != 1. || optimal_ϕ_function(ϕ_prop) < 0
            ϕ_n = fzero(optimal_ϕ_function, [ϕ_n1, ϕ_prop], xtol = 0.)
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
    # if !tempered_update
        incremental_weights = exp.((ϕ_n - ϕ_n1)*get_loglh(cloud))
    # else
        # incremental_weights = exp.((ϕ_n1 - ϕ_n)*revised_loglh + (ϕ_n - ϕ_n1)*new_loglh)
    # end

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
    push!(cloud.ESS, 1/sum(normalized_weights.^2))

    # If this assertion does not hold then there are likely too few particles
    @assert !isnan(cloud.ESS[i]) "no particles have non-zero weight"

    # Resample if the degeneracy/effective sample size metric falls below the accepted threshold
    if (cloud.ESS[i] < threshold)
        new_inds = resample(normalized_weights; method = get_setting(m, :resampler_smc),
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

    # Calculate the adaptive c-step to be used as a scaling coefficient in the mutation MH step
    c = c*(0.95 + 0.10*exp(16*(cloud.accept - target))/(1 + exp(16*(cloud.accept - target))))
    m <= Setting(:step_size_smc, c)
    cloud.c = c

    if parallel
        R = cov(cloud)
        new_particles = @parallel (vcat) for j in 1:n_parts
            mutation(m, data, cloud.particles[j], R, ϕ_n, ϕ_n1; c = c)
        end
    else
        R = cov(cloud)
        new_particles = [mutation(m, data, cloud.particles[j], R, ϕ_n, ϕ_n1; c = c) for j = 1:n_parts]
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

    if !m.testing && !tempered_update
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
    elseif tempered_update
        new_vintage = get_setting(m, :updated_data_vintage)
        jldopen(rawpath(m, "estimate", "smc_cloud_updvint=$new_vintage.jld"), "r+") do file
            write(file, "up_cloud", cloud)
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

# For doing a combination of data tempering and likelihood tempering to incorporate
# new data into particle cloud sample
# Cloud is a obtained from loading from a .jld an old particle cloud pertaining to the old data matrix
# from a previous run of SMC, located in the relevant path.
# new_data is a new data matrix corresponding to both revised and newly updated periods of data to be
# incorporated into creating a new cloud/sample.
function smc{S<:AbstractFloat}(m::AbstractModel, cloud::ParticleCloud, old_data::Matrix{S},
                               new_data::Matrix{S}; verbose::Symbol=:low)
    new_vintage = get_setting(m, :updated_data_vintage)
    original_data_range  = quarter_range(date_mainsample_start(m), date_mainsample_end(m))
    updated_data_range   = quarter_range(get_setting(m, :date_updatedsample_start), get_setting(m, :date_updatedsample_end))
    revised_quarters = intersect(original_data_range, updated_data_range)
    new_quarters = setdiff(updated_data_range, original_data_range)

    n_particles = length(cloud)
    unchanged_loglh = zeros(n_particles)
    revised_loglh   = zeros(n_particles)
    new_loglh       = zeros(n_particles)

    for (i,p) in enumerate(cloud.particles)
        update!(m, p.value)
        unchanged_loglh[i] = likelihood(m, old_data[:, 1:end-length(revised_quarters)])
        revised_loglh[i]   = likelihood(m, old_data[:, end-length(revised_quarters)+1:end])
        new_loglh[i]       = likelihood(m, new_data)
    end

    data = hcat(old_data, new_data)

    jldopen(rawpath(m, "estimate", "smc_cloud_updvint=$new_vintage.jld"), "w") do file
        write(file, "old_cloud", cloud)
        write(file, "unchanged_loglh", unchanged_loglh)
        write(file, "revised_loglh", revised_loglh)
        write(file, "new_loglh", new_loglh)
    end

    smc(m, data; verbose=verbose, tempered_update=true)
end
