"""
```
smc(m::AbstractModel, data::Matrix)
smc(m::AbstractModel, data::DataFrame)
smc(m::AbstractModel)
```

### Arguments:

- `m`: A model object, from which its parameters values, prior dists, and bounds will be referenced
- `data`: A matrix or data frame containing time series of the observables to be used in the calculation of the posterior/likelihood

### Keyword Arguments:
- `verbose`: The desired frequency of function progress messages printed to standard out.
	- `:none`: No status updates will be reported.
	- `:low`: Status updates for SMC initialization and recursion will be included.
	    - `:high`: Status updates for every iteration of SMC is output, which includes the mu and sigma of each individual parameter draw after each iteration, as well as the calculated acceptance rate, the ESS, and the number of times resampled.
- `start_approx_post`: Have the initial parameters drawn from the approximated posterior (a degenerate multivariate normal
centered at the mode).

### Outputs

- `mu`: The mean of a particular parameter as calculated by taking the weighted sum of that parameter from the particle draws made in that particular stage Nφ.
- `sig`: The standard deviation of a particular parameter.

### Overview

Sequential Monte Carlo can be used in lieu of Random Walk Metropolis Hastings to generate parameter samples from high-dimensional parameter spaces using sequentially constructed proposal densities.

The implementation here is based upon Edward Herbst and Frank Schorfheide's 2014 paper 'Sequential Monte Carlo Sampling for DSGE Models' and the code accompanying their book 'Bayesian Estimation of DSGE Models'.

SMC is broken up into three main steps:

- `Correction`: Reweight the particles from stage n-1 by defining "incremental weights", incweight, which gradually incorporate the likelihood function p(Y|θ(i, n-1)) into the particle weights.
- `Selection`: Resample the particles if the distribution of particles begins to degenerate, according to a tolerance level ESS < N/2. The current resampling technique employed is systematic resampling.

- `Mutation`: Propagate particles {θ(i), W(n)} via N(MH) steps of a Metropolis Hastings algorithm.
"""
function smc(m::AbstractModel, data::Matrix; verbose::Symbol = :low, new_data::Matrix = Array{Float64,2}())
    ########################################################################################
    ### Setting Parameters
    ########################################################################################

    parallel = get_setting(m, :use_parallel_workers)

    n_parts = get_setting(m, :n_particles)
    n_params = n_parameters(m)
    n_blocks = get_setting(m, :n_smc_blocks)
    n_Φ = get_setting(m, :n_Φ)

    c = get_setting(m, :step_size_smc)
    accept = get_setting(m, :init_accept)
    target = get_setting(m, :target_accept)
    λ = get_setting(m, :λ)

    # Instantiating ParticleCloud object
    cloud = ParticleCloud(m, n_parts)
    cloud.n_Φ = n_Φ
    cloud.c = c
    cloud.tempering_schedule = ((collect(1:n_Φ)-1)/(n_Φ-1)).^λ
    cloud.accept = accept

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

    # Particle draws from the parameter's marginal priors
    # Modifies the cloud object in place to update draws, loglh, & logpost
    if isempty(new_data)
        initial_draw(m, data, cloud)
    else
        cloud = 
    end

    if m.testing
        writecsv("draws.csv", get_vals(cloud))
    end

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        init_stage_print(cloud;verbose=verbose)
    end

    ########################################################################################
    ### Recursion
    ########################################################################################

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("\n\n SMC Recursion starts \n\n")
    end

    # Draws saved for debugging purposes
    # Saving all clouds for debugging purposes
    # simfile = jldopen(rawpath(m, "estimate", "smc_cloud.jld"), "w")
    # cloud_store = g_create(simfile, "cloud_store")
    # cloud_store["cloud_1"] = cloud

    total_sampling_time = 0.

    for i = 2:n_Φ
        tic()
        cloud.stage_index = i

    ########################################################################################
    ### Step 1: Correction
    ########################################################################################

	# Incremental weights
    if isempty(new_data)
        inc_weight = exp((cloud.tempering_schedule[i] - cloud.tempering_schedule[i-1]) * get_loglh(cloud))
    else
        inc_weight = 
    end

    # Update weights
    update_weights!(cloud, inc_weight)

	# Normalize weights
    normalize_weights!(cloud)

    if m.testing
        writecsv("weights.csv", get_weights(cloud))
    end

    ########################################################################################
    ### Step 2: Selection
    ########################################################################################

    cloud.ESS = 1/sum(get_weights(cloud).^2)
    @assert !isnan(cloud.ESS) "no particles have non-zero weight"
    if (cloud.ESS < n_parts/2)
        new_inds = resample(get_weights(cloud); method = get_setting(m, :resampler_smc),
                            parallel = parallel, testing = m.testing)
        # update parameters/logpost/loglh with resampled values
        # reset the weights to 1/n_parts
        cloud.particles = [deepcopy(cloud.particles[i]) for i in new_inds]
        reset_weights!(cloud)
        cloud.resamples += 1
    end

    ########################################################################################
    ### Step 3: Mutation
    ########################################################################################

    c = c*(0.95 + 0.10*exp(16*(cloud.accept - target))/(1 + exp(16*(cloud.accept - target))))
    # Updating c step size in model object and cloud object
    m <= Setting(:step_size_smc, c)
    cloud.c = c

    if parallel
        R = cov(cloud)
        current_φ = cloud.tempering_schedule[i]
        previous_φ = cloud.tempering_schedule[i-1]
        #out = pmap(j -> mutation_block_RWMH(m, data, cloud.particles[j], R, current_φ, previous_φ, 1:n_parts))
        new_particles = @parallel (vcat) for j in 1:n_parts
            mutation(m, data, cloud.particles[j], R, current_φ, previous_φ)
        end
    else
        R = cov(cloud)
        current_φ = cloud.tempering_schedule[i]
        previous_φ = cloud.tempering_schedule[i-1]
        new_particles = [mutation(m, data, cloud.particles[j], R, current_φ, previous_φ)  for j = 1:n_parts]
    end

    cloud.particles = new_particles
    update_acceptance_rate!(cloud) # update average acceptance rate

    ########################################################################################
    ### Timekeeping and Output Generation
    ########################################################################################

    total_sampling_time += toq()

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        end_stage_print(cloud, total_sampling_time; verbose=verbose)
    end

    ########################################################################################
    ### Saving data
    ########################################################################################

    # Draws saved for debugging purposes
    # Save both the standard format .h5 file of parameters by particles for reading
    # into the forecast step, but also save the last cloud object as a .jld
    # cloud_store["cloud_$i"] = cloud

    if m.testing
        writecsv("mutated_values.csv",get_vals(cloud))
        break
    end

    end

    if !m.testing
        simfile = h5open(rawpath(m,"estimate","smcsave.h5"),"w")
        particle_store = d_create(simfile, "smcparams", datatype(Float32),
                                  dataspace(n_parts, n_params))
        for i in 1:length(cloud)
            particle_store[i,:] = cloud.particles[i].value
        end
        close(simfile)
        save(rawpath(m,"estimate","smc_cloud.jld"),"cloud",cloud)
    end
    # Draws saved for debugging purposes
    # close(simfile)
    # Access saved clouds by jldopen then read(file,"cloud_store")
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


"""
```
initial_draw(m::AbstractModel, data::Array{Float64}, c::ParticleCloud)
```

Draw from a general starting distribution (set by default to be from the prior) to initialize the SMC algorithm.
Returns a tuple (logpost, loglh) and modifies the particle objects in the particle cloud in place.

"""
function initial_draw(m::AbstractModel, data::Array{Float64}, c::ParticleCloud)
    dist_type = get_setting(m, :initial_draw_source)
    if dist_type == :normal
        params = zeros(n_parameters(m))
        hessian = zeros(n_parameters(m),n_parameters(m))
        try
            file = h5open(rawpath(m, "estimate", "paramsmode.h5"), "r")
            params = read(file,"params")
            close(file)
        catch
            throw("There does not exist a valid mode file at "*rawpath(m,"estimate","paramsmode.h5"))
        end

        try
            file = h5open(rawpath(m, "estimate", "hessian.h5"), "r")
            hessian = read(file, "hessian")
            close(file)
        catch
            throw("There does not exist a valid hessian file at "*rawpath(m,"estimate","hessian.h5"))
        end

        S_diag, U = eig(hessian)
        big_eig_vals = find(x -> x>1e-6, S_diag)
        rank = length(big_eig_vals)
        n = length(params)
        S_inv = zeros(n, n)
        for i = (n-rank+1):n
            S_inv[i, i] = 1/S_diag[i]
        end
        hessian_inv = U*sqrt(S_inv)
        dist = DSGE.DegenerateMvNormal(params, hessian_inv)
    end

    n_part = length(c)
    draws =
    dist_type == :prior ? rand(m.parameters, n_part) : rand(dist, n_part)

    loglh = zeros(n_part)
    logpost = zeros(n_part)
    for i in 1:size(draws)[2]
        success = false
        while !success
            try
                update!(m, draws[:, i])
                loglh[i] = likelihood(m, data)
                logpost[i] = prior(m)
            catch
                draws[:, i] =
                dist_type == :prior ? rand(m.parameters, 1) : rand(dist, 1)
                continue
            end
            success = true
        end
    end
    update_draws!(c, draws)
    update_loglh!(c, loglh)
    update_logpost!(c, logpost)
end

function init_stage_print(cloud::ParticleCloud; verbose::Symbol=:low)
    μ = weighted_mean(cloud)
    σ = weighted_std(cloud)
	println("--------------------------")
        println("Iteration = $(cloud.stage_index) / $(cloud.n_Φ)")
	println("--------------------------")
        println("phi = $(cloud.tempering_schedule[cloud.stage_index])")
	println("--------------------------")
        println("c = $(cloud.c)")
	println("--------------------------")
    if VERBOSITY[verbose] >= VERBOSITY[:high]
        for n=1:length(cloud.particles[1])
            println("$(cloud.particles[1].keys[n]) = $(round(μ[n], 5)), $(round(σ[n], 5))")
	    end
    end
end

function end_stage_print(cloud::ParticleCloud, total_sampling_time::Float64; verbose::Symbol=:low)
    total_sampling_time_minutes = total_sampling_time/60
    expected_time_remaining_sec = (total_sampling_time/cloud.stage_index)*(cloud.n_Φ - cloud.stage_index)
    expected_time_remaining_minutes = expected_time_remaining_sec/60

    μ = weighted_mean(cloud)
    σ = weighted_std(cloud)
    println("--------------------------")
        println("Iteration = $(cloud.stage_index) / $(cloud.n_Φ)")
        println("time elapsed: $(round(total_sampling_time_minutes, 4)) minutes")
        println("estimated time remaining: $(round(expected_time_remaining_minutes, 4)) minutes")
    println("--------------------------")
        println("phi = $(cloud.tempering_schedule[cloud.stage_index])")
    println("--------------------------")
        println("c = $(cloud.c)")
        println("accept = $(cloud.accept)")
        println("ESS = $(cloud.ESS)   ($(cloud.resamples) total resamples.)")
    println("--------------------------")
    if VERBOSITY[verbose] >= VERBOSITY[:high]
        for n=1:length(cloud.particles[1])
            println("$(cloud.particles[1].keys[n]) = $(round(μ[n], 5)), $(round(σ[n], 5))")
        end
    end
end
