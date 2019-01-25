"""
```
smc(m::AbstractModel, data::Matrix; verbose::Symbol, old_data::Matrix)
smc(m::AbstractModel, data::DataFrame)
smc(m::AbstractModel)
```

### Arguments:

- `m`: A model object, which stores parameter values, prior dists, bounds, and various other settings that will be referenced
- `data`: A matrix or dataframe containing the time series of the observables used in the calculation of the posterior/likelihood
- `old_data`: A matrix containing the time series of observables of previous data (with `data` being the new data) for the purposes of a time tempered estimation (that is, using the posterior draws from a previous estimation as the initial set of draws for an estimation with new data)

### Keyword Arguments:
- `verbose`: The desired frequency of function progress messages printed to standard out.
	- `:none`: No status updates will be reported.
	- `:low`: Status updates for SMC initialization and recursion will be included.
	- `:high`: Status updates for every iteration of SMC is output, which includes the mean and standard deviation of each individual parameter draw after each iteration, as well as the calculated acceptance rate, the ESS, and the number of times resampled.

### Outputs

- `cloud`: The ParticleCloud object containing all of the information about the parameter values from the sample, their respective log-likelihoods, the ESS schedule, tempering schedule etc., which is saved in the saveroot.

### Overview

Sequential Monte Carlo can be used in lieu of Random Walk Metropolis Hastings to generate parameter samples from high-dimensional parameter spaces using sequentially constructed proposal densities to be used in iterative importance sampling.

The implementation here is based on Edward Herbst and Frank Schorfheide's 2014 paper 'Sequential Monte Carlo Sampling for DSGE Models' and the code accompanying their book 'Bayesian Estimation of DSGE Models'.

SMC is broken up into three main steps:

- `Correction`: Reweight the particles from stage n-1 by defining incremental weights, which gradually "temper in" the likelihood function p(Y|θ)^(ϕ_n - ϕ_n-1) into the normalized particle weights.
- `Selection`: Resample the particles if the distribution of particles begins to degenerate, according to a tolerance level for the ESS.
- `Mutation`: Propagate particles {θ(i), W(n)} via N(MH) steps of a Metropolis Hastings algorithm.
"""
function smc(m::AbstractModel, data::Matrix{Float64};
             verbose::Symbol = :low, old_data::Matrix{Float64} = Matrix{Float64}(size(data, 1), 0),
             old_cloud::ParticleCloud = ParticleCloud(m, 0))
    ########################################################################################
    ### Setting Parameters
    ########################################################################################

    # RECA
    fort_file    = "smc-sw-new-mix-npart-12000-nintmh-1-nphi-500-prior-b1-trial1-phibend-jan24-mixron/"
    fortran_path = "/home/rcerxs30/SLICOT-2018-12-19/dsge-smc/fortran/" * fort_file
    #"/data/dsge_data_dir/dsgejl/reca/SMCProject/specfiles/fortran/"
    #fortESS      = vec(readdlm(fortran_path * "ESS.txt"))
    stepprobs    = vec(readdlm(fortran_path * "stepprobstemp.txt"))

    # General
    parallel = get_setting(m, :use_parallel_workers)
    n_parts  = get_setting(m, :n_particles)
    n_params = n_parameters(m)
    n_blocks = get_setting(m, :n_smc_blocks)
    n_steps  = get_setting(m, :n_mh_steps_smc)

    use_chand_recursion = get_setting(m, :use_chand_recursion)

    if any(isnan.(data)) & use_chand_recursion
        error("Cannot use Chandrasekhar recursions with missing data")
    end

    # Time Tempering
    tempered_update = !isempty(old_data)

    # Quick check that if there is a tempered update that the old vintage and current vintage are different
    if tempered_update
        old_vintage = get_setting(m, :previous_data_vintage)
        @assert old_vintage != data_vintage(m)
    end

    # Step 0 (ϕ schedule) settings
    i = 1                           # The index tracking the stage of the algorithm
    j = 2                           # The index tracking to the fixed_schedule entry that ϕ_prop is set as
    resampled_last_period = false   # To ensure proper resetting of ESS_bar right after resample
    ϕ_n = 0.                        # Instantiating ϕ_n and ϕ_prop variables to be referenced in their
    ϕ_prop = 0.                     # respective while loop conditions
    use_fixed_schedule = get_setting(m, :adaptive_tempering_target_smc) == 0.0 #get_setting(m, :use_fixed_schedule)
    λ = get_setting(m, :λ)
    n_Φ = get_setting(m, :n_Φ)
    tempering_target = get_setting(m, :adaptive_tempering_target_smc)

    # Step 2 (Correction) settings
    resampling_method = get_setting(m, :resampler_smc)
    threshold_ratio   = get_setting(m, :resampling_threshold)
    threshold         = threshold_ratio * n_parts

    # Step 3 (Mutation) settings
    c      = get_setting(m, :step_size_smc)
    target = accept = get_setting(m, :target_accept)
    α      = get_setting(m, :mixture_proportion)

    fixed_para_inds = find([θ.fixed for θ in m.parameters])
    free_para_inds  = find([!θ.fixed for θ in m.parameters])
    n_free_para     = length(free_para_inds)

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

    if tempered_update
        if isempty(old_cloud)
            # Load the previous ParticleCloud as the starting point for time tempering
            loadpath = rawpath(m, "estimate", "smc_cloud.jld")
            #loadpath = rawpath(m, "estimate", "smc_cloud.jld", ["adpt="*string(tempering_target)])
            loadpath = replace(loadpath, r"vint=[0-9]{6}", "vint="*old_vintage)

            cloud = load(loadpath, "cloud")
        else
            cloud = old_cloud
        end
        initialize_cloud_settings!(m, cloud; tempered_update = tempered_update)
        initialize_likelihoods!(m, data, cloud, parallel = parallel, verbose = verbose)
    else
        # Instantiating ParticleCloud object
        cloud = ParticleCloud(m, n_parts)

        # Modifies the cloud object in place to update draws, loglh, & logpost
        initial_draw!(m, data, cloud, fortran_path, parallel = parallel,
                      use_chand_recursion = use_chand_recursion,
                      verbose = verbose) # RECA edited this

        initialize_cloud_settings!(m, cloud; tempered_update = tempered_update)
    end

    # Fixed schedule for construction of ϕ_prop
    if use_fixed_schedule
        cloud.tempering_schedule = ((collect(1:n_Φ)-1)/(n_Φ-1)).^λ
    else
        proposed_fixed_schedule = ((collect(1:n_Φ)-1)/(n_Φ-1)).^λ
    end

    # Instantiate incremental and normalized weight matrices to be used for logMDD calculation
    w_matrix = zeros(n_parts, 1)
    if tempered_update
        W_matrix = similar(w_matrix)
        for k in 1:n_parts
            W_matrix[k] = cloud.particles[k].weight
        end
    else
        W_matrix = fill(1/n_parts, (n_parts,1))
    end
    z_matrix = ones(1)

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        init_stage_print(cloud; verbose = verbose, use_fixed_schedule = use_fixed_schedule)
    end

    ########################################################################################
    ### Recursion
    ########################################################################################
    # RECA
    resample_periods = get_resample_periods(fortran_path)

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("\n\n SMC recursion starts \n\n")
    end

    while ϕ_n < 1.

        tic()
        cloud.stage_index = i += 1

        #RECA
        para_cur = readdlm(fortran_path * convert_string(i-1) * "parasim.txt")
        @test_matrix_approx_eq get_vals(cloud) transpose(para_cur)

        if !isapprox(get_vals(cloud), transpose(para_cur))
            @show maximum(abs.(get_vals(cloud) - transpose(para_cur)) ./ transpose(para_cur))
        end
        #END RECA

        ########################################################################################
        ### Step 0: Setting ϕ_n (either adaptively or by the fixed schedule)
        ########################################################################################
        ϕ_n1 = cloud.tempering_schedule[i-1]

        if use_fixed_schedule
            ϕ_n = cloud.tempering_schedule[i]
        else
            ϕ_n, resampled_last_period, j, ϕ_prop = solve_adaptive_ϕ(cloud, proposed_fixed_schedule,
                                                                     i, j, ϕ_prop, ϕ_n1,
                                                                     tempering_target,
                                                                     resampled_last_period)
        end

        ########################################################################################
        ### Step 1: Correction
        ########################################################################################

        # Calculate incremental weights (if no old data, get_old_loglh(cloud) returns zero)
        incremental_weights = exp.((ϕ_n1 - ϕ_n)*get_old_loglh(cloud) + (ϕ_n - ϕ_n1)*get_loglh(cloud))

        # Update weights
        update_weights!(cloud, incremental_weights)
        mult_weights = get_weights(cloud)

        # Normalize weights
        normalize_weights!(cloud)
        normalized_weights = get_weights(cloud)

        # RECA
        fort_wtsim = vec(readdlm(fortran_path * convert_string(i) * "cor_weights.txt"))
        @assert isapprox(normalized_weights, fort_wtsim, nans=true, atol=1e-6)

        push!(z_matrix, sum(mult_weights))
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
            if i in resample_periods # RECA: test to ensure resample at same periods
                print_with_color(:blue, "Resample occurs in FORTRAN set of resample periods.\n")
                new_inds = resample(i, fortran_path)

                # update parameters/logpost/loglh with resampled values
                # reset the weights to 1/n_parts
                cloud.particles = [deepcopy(cloud.particles[i]) for i in new_inds]
                reset_weights!(cloud)
                cloud.resamples += 1
                resampled_last_period = true
                W_matrix[:, i] = fill(1/n_parts, (n_parts,1))
            else
                print_with_color(:red, "RESAMPLE NOT IN FORTRAN SET OF RESAMPLE PERIODS - ESS TOO LOW.\n")
                #new_inds = resample(normalized_weights; method = resampling_method)
            end
        end

        ########################################################################################
        ### Step 3: Mutation
        ########################################################################################

        # Calculate the adaptive c-step to be used as a scaling coefficient in the mutation MH step
        c = c*(0.95 + 0.10*exp(16*(cloud.accept - target))/(1. + exp(16.*(cloud.accept - target))))
        cloud.c = c

        θ_bar = weighted_mean(cloud)
        R     = weighted_cov(cloud)

        # RECA: Testing for validity of run
        mu       = vec(readdlm(fortran_path * convert_string(i) * "mean.txt"))
        fort_var = readdlm(fortran_path * convert_string(i) * "var.txt")
        fort_c   = readdlm(fortran_path * convert_string(i) * "scale.txt")[1]

        if !(c ≈ fort_c) @show c, fort_c end
        @assert c ≈ fort_c
        if !(θ_bar ≈ mu)
            print_with_color(:red, "θ_bar FAILS TO BE EQUAL TO mu WITHIN 1e-6.\n")
            @show θ_bar, mu
        end
        @assert isapprox(θ_bar, mu, atol=1e-3) #θ_bar ≈ mu

        #@show cloud.ESS[i] ≈ fortESS[i], cloud.ESS[i], fortESS[i]

        # add to itself and divide by 2 to ensure marix is positive semi-definite symmetric
        # (not off due to numerical error) and values haven't changed
        R_fr = (R[free_para_inds, free_para_inds] + R[free_para_inds, free_para_inds]') / 2

        # Confirm these babies are the same
        @test_matrix_approx_eq R_fr fort_var[free_para_inds, free_para_inds]
        if !isapprox(R_fr, fort_var[free_para_inds, free_para_inds])
          @show sum(abs.(R_fr - fort_var[free_para_inds, free_para_inds])./fort_var[free_para_inds, free_para_inds])
          @show maximum(abs.(R_fr-fort_var[free_para_inds,free_para_inds])./fort_var[free_para_inds,free_para_inds])
        end

        # MvNormal centered at ̄θ with var-cov ̄Σ, subsetting out the fixed parameters
        d = MvNormal(θ_bar[free_para_inds], R_fr)

        # New way of generating blocks
        # blocks_free = generate_free_blocks(n_free_para, n_blocks)
        # blocks_all  = generate_all_blocks(blocks_free, free_para_inds)
        blocks_all   = generate_all_blocks(n_free_para, n_blocks, fortran_path, i) # RECA: fn names
        blocks_free  = generate_free_blocks(blocks_all, free_para_inds) # RECA

        # RECA: mixr gives which distribution we ought draw from
        mixr     = readdlm(fortran_path * convert_string(i) * "mixr.txt")
        eps      = readdlm(fortran_path * convert_string(i) * "eps.txt")
        step_pr  = readdlm(fortran_path * convert_string(i) * "step_pr.txt")
        step_lik = readdlm(fortran_path * convert_string(i) * "step_lik.txt")
        step_p1  = readdlm(fortran_path * convert_string(i) * "step_p1.txt")
        alp      = readdlm(fortran_path * convert_string(i) * "step_alp.txt")
        fortpara = readdlm(fortran_path * convert_string(i) * "parasim.txt")
        fortpost = readdlm(fortran_path * convert_string(i) * "postsim.txt")
        fortlik  = readdlm(fortran_path * convert_string(i) * "liksim.txt")
        bvar     = readdlm(fortran_path * convert_string(i) * "bvar.txt")


        if parallel
            new_particles = @parallel (vcat) for j in 1:n_parts
                # RECA: Testing against FORTRAN
                # code below returns the same ESS value, very comparable logMDD
                mutation(m, data, cloud.particles[j], d, blocks_free, blocks_all, ϕ_n, ϕ_n1;
                         c = c, α = α, old_data = old_data,
                         use_chand_recursion = use_chand_recursion, verbose = verbose,
                         # RECA'S TESTING VARIABLES:
                         mixr = mixr[j,:],
                         eps = eps[(j-1) * n_steps + 1:j * n_steps, :],
                         stepprobs = stepprobs[(i-2) * n_parts * n_steps * n_blocks + (j-1) * n_steps * n_blocks + 1:(i-2) * n_parts * n_steps * n_blocks + (j-1) * n_steps * n_blocks + n_blocks * n_steps],
                         step_pr = step_pr[j,:], step_lik = step_lik[j,:],
                         step_p1 = step_p1[j,:], alp = alp[j,:], mu = mu, bvar = bvar)
            end
        else
            new_particles = [mutation(m, data, cloud.particles[j], d, blocks_free, blocks_all,
                                      ϕ_n, ϕ_n1; c = c, α = α, old_data = old_data,
                                      use_chand_recursion = use_chand_recursion, verbose = verbose,
                                      # RECA'S TESTING VARIABLES:
                                      pnum = j,
                                      mixr = mixr[j,:],
                                      eps = eps[(j-1) * n_steps + 1:j * n_steps, :],
                                      stepprobs = stepprobs[(i-2) * n_parts * n_steps * n_blocks + (j-1) * n_steps * n_blocks + 1:(i-2) * n_parts * n_steps * n_blocks + (j-1) * n_steps * n_blocks + n_blocks * n_steps],
                                      step_pr = step_pr[j,:], step_lik = step_lik[j,:],
                                      step_p1 = step_p1[j,:], alp = alp[j,:],
                                      mu = mu, bvar = bvar) for j=1:n_parts]
        end

        cloud.particles = new_particles
        update_acceptance_rate!(cloud) # update average acceptance rate

        ########################################################################################
        ### Reca: Testing Against FORTRAN
        ########################################################################################
        for kk=1:length(vec(fortlik))
            if !((abs(get_loglh(cloud)[kk] - vec(fortlik)[kk]) / abs(vec(fortlik)[kk])) < 1e-3)
                @show kk, abs(get_loglh(cloud)[kk] - vec(fortlik)[kk]) / abs(vec(fortlik)[kk]), vec(fortlik)[kk], get_loglh(cloud)[kk]
            end
#            sleep(0.5)
            @assert (abs(get_loglh(cloud)[kk] - vec(fortlik)[kk]) / abs(vec(fortlik)[kk])) < 1e-2
        end
        @test_matrix_approx_eq get_vals(cloud) transpose(fortpara)
        if !isapprox(get_vals(cloud), transpose(fortpara))
            @show !isapprox(get_vals(cloud), transpose(fortpara), atol=1e-3)
            @show maximum(abs.(get_vals(cloud) - transpose(fortpara)) ./ transpose(fortpara))
        end
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
        simfile = h5open(rawpath(m, "estimate", "smcsave.h5"), "w")
        #simfile = h5open(rawpath(m, "estimate", "smcsave.h5", ["adpt="*string(tempering_target)]),"w")
        particle_store = d_create(simfile, "smcparams", datatype(Float32),
                                  dataspace(n_parts, n_params))
        for i in 1:length(cloud)
            particle_store[i,:] = cloud.particles[i].value
        end
        close(simfile)
        #jldopen(rawpath(m, "estimate", "smc_cloud.jld", ["adpt="*string(tempering_target)]), "w") do file
        jldopen(rawpath(m, "estimate", "smc_cloud.jld"), "w") do file
            write(file, "cloud", cloud)
            write(file, "w", w_matrix)
            write(file, "W", W_matrix)
            write(file, "z", z_matrix)
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
