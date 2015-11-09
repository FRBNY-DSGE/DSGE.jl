# This program produces and saves draws from the posterior distribution of
# the parameters.

#=
doc """
estimate{T<:AbstractDSGEModel}(m::T; verbose::Symbol=:low, proposal_covariance=[])

### Parameters:
- `m`: The model object

### Optional Arguments:
- `verbose`: The desired frequency of function progress messages printed to standard out.

   - `:none`: No status updates will be reported.

   - `:low`: Status updates will be provided in csminwel and at each block in Metropolis-Hastings.

   - `:high`: Status updates provided at each iteration in Metropolis-Hastings.

- `proposal_covariance`: Used to test the metropolis_hastings algorithm with a precomputed covariance matrix for the proposal distribution. When the Hessian is singular, eigenvectors corresponding to zero eigenvectors are not well defined, so eigenvalue decomposition can cause problems. Passing a precomputed matrix allows us to ensure that the rest of the routine has not broken.

### Description
This routine implements the full estimation stage of the FRBNY DSGE model.
"""->
=#
function estimate{T<:AbstractDSGEModel}(m::T; verbose::Symbol=:low, proposal_covariance=[])

    ########################################################################################
    ### Step 1: Initialize
    ########################################################################################

    # Load data
    h5 = h5open(inpath(m, "data","data_$(data_vintage(m)).h5"), "r")
    YY = read(h5["YY"])
    close(h5)

    post = posterior(m, YY)[:post]


    ########################################################################################
    ### Step 2: Find posterior mode (if reoptimizing, run csminwel)
    ########################################################################################
    
    # Specify starting mode

    if VERBOSITY[verbose] > VERBOSITY[:none] 
        println("Reading in previous mode")
    end
    
    mode = []

    if reoptimize(m)
        h5 = h5open(inpath(m, "user", "mode_in.h5"),"r")
        mode = read(h5["params"])   #it's mode in mode_in_optimized, but params in mode_in
        close(h5)
    else
        h5 = h5open(inpath(m, "user", "mode_in_optimized.h5"),"r")
        mode = read(h5["mode"])
        close(h5)
    end


    update!(m, mode)

    if reoptimize(m)
        println("Reoptimizing...")
        
        # Inputs to minimization algorithm
        function posterior_min!{T<:AbstractFloat}(x::Vector{T})
            tomodel!(m,x)
            return -posterior(m, YY; catch_errors=true)[:post]
        end

        xh = toreal(m.parameters)
        H = 1e-4 * eye(num_parameters(m))
        nit = 1000
        crit = 1e-10
        converged = false

        # If the algorithm stops only because we have exceeded the maximum number of
        # iterations, continue improving guess of modal parameters
        while !converged
            out, H = csminwel(posterior_min!, xh, H; model=m, ftol=crit, iterations=nit, show_trace=true, verbose=verbose)
            xh = out.minimum
            converged = !out.iteration_converged
        end

        # Transform modal parameters so they are no longer bounded (i.e., allowed
        # to lie anywhere on the real line).
        tomodel!(m, xh)
        mode = [param.value for param in m.parameters]

        # Write mode to file
        h5 = h5open(rawpath(m, "estimate", "mode_out.h5"),"w")
        h5["mode"] = mode
        close(h5)
        
    end

    
    ########################################################################################
    ### Step 3: Compute proposal distribution
    ###
    ### In Metropolis-Hastings, we draw sample parameter vectors from
    ### the proposal distribution, which is a degenerate multivariate
    ### normal centered at the mode. Its variance is the inverse of
    ### the hessian. We find the inverse via eigenvalue decomposition.
    ########################################################################################

    # Calculate the Hessian at the posterior mode
    hessian = if recalculate_hessian(m)
        if VERBOSITY[verbose] > VERBOSITY[:none] 
            println("Recalculating Hessian...")
        end
        
        hessian, _ = hessian!(m, mode, YY; verbose=true)

        h5open(rawpath(m, "estimate","hessian.h5"),"w") do file
            file["hessian"] = hessian
        end

        hessian

    # Read in a pre-optimized mode
    else
        if VERBOSITY[verbose] > VERBOSITY[:none]
            println("Using pre-calculated Hessian")
        end

        hessian = h5open(inpath(m, "user", "hessian_optimized.h5"),"r") do file
            read(file, "hessian")
        end
    end

    # Compute inverse hessian and create proposal distribution, or
    # just create it with the given cov matrix if we have it
    propdist = if length(proposal_covariance) == 0
        # Make sure the mode and hessian have the same number of parameters
        n = length(mode)
        @assert (n, n) == size(hessian)

        # Compute the inverse of the Hessian via eigenvalue decomposition
        S_diag, U = eig(hessian)
        big_eigvals = find(x -> x > 1e-6, S_diag)
        rank = length(big_eigvals)
    
        S_inv = zeros(n, n)

        for i = (n-rank+1):n
            S_inv[i, i] = 1/S_diag[i]
        end
        
        hessian_inv = U*sqrt(S_inv) #this is the inverse of the hessian
        DegenerateMvNormal(mode, hessian_inv, rank)
    else
        DegenerateMvNormal(mode, proposal_covariance)
    end
    
    if propdist.rank != num_parameters_free(m)
        println("problem –    shutting down dimensions")
    end


    ########################################################################################
    ### Step 4: Sample from posterior using Metropolis-Hastings algorithm
    ########################################################################################
    
    # Set the jump size for sampling
    cc0 = 0.01
    cc = 0.09

    metropolis_hastings(propdist, m, YY, cc0, cc; verbose=verbose);

    
    ########################################################################################
    ### Step 5: Calculate and save parameter covariance matrix
    ########################################################################################

    compute_parameter_covariance(m);

    return nothing 
end


#=
doc"""
metropolis_hastings{T<:AbstractFloat}(propdist::Distribution, m::AbstractDSGEModel, YY::Matrix{T}, cc0::T, cc::T; verbose::Symbol = :low, fix_seed::Bool=true)

### Parameters
* `propdist` The proposal distribution that Metropolis-Hastings begins sampling from.
* `m`: The model object
* `YY`: Data matrix for observables
* `cc0`: Jump size for initializing Metropolis-Hastings.
* `cc`: Jump size for the rest of Metropolis-Hastings.

### Optional Arguments
* `verbose`: The desired frequency of function progress messages printed to standard out.

   - `:none`: No status updates will be reported.

   - `:low`: Status updates provided at each block.

   - `:high`: Status updates provided at each draw.

* `fix_seed`: fix the seed of the random number generator

### Description
Implements the Metropolis-Hastings MCMC algorithm for sampling from the posterior distribution of the parameters.
"""
=#
function metropolis_hastings{T<:AbstractFloat}(propdist::Distribution, m::AbstractDSGEModel,
       YY::Matrix{T}, cc0::T, cc::T; verbose::Symbol=:low, fix_seed::Bool=true)

    # If testing, set the random seeds at fixed numbers
    
    if fix_seed || m.testing
        srand(m.rng, 654)
    end
        
    # Set number of draws, how many we will save, and how many we will burn
    # (initialized here for scoping; will re-initialize in the while loop)

    n_blocks = 0
    n_sim = 0
    n_times = 0
    n_burn = 0
    n_params = num_parameters(m)

    # Initialize algorithm by drawing para_old from a normal distribution centered on the
    # posterior mode until the parameters are within bounds or the posterior value is sufficiently large.
    para_old = rand(propdist, m; cc=cc0)
    post_old = -Inf
    like_old = -Inf

    TTT_old = []
    RRR_old = []
    CCC_old = []

    zend_old = []
    ZZ_old = []
    DD_old = []
    QQ_old = []

    initialized = false

    while !initialized

        n_blocks = num_mh_blocks(m)
        n_sim    = num_mh_simulations(m)
        n_burn   = num_mh_burn(m)
        n_times  = mh_thinning_step(m)

        post_out = posterior!(m, para_old, YY; mh=true)
        post_old, like_old, out = post_out[:post], post_out[:like], post_out[:mats]
        
        if post_old > -Inf
            propdist.μ = para_old

            TTT_old = out[:TTT]
            RRR_old = out[:RRR]
            CCC_old = out[:CCC]
            zend_old  = out[:zend]

            initialized = true
        end

    end

    # Report number of blocks that will be used 
    if VERBOSITY[verbose] > VERBOSITY[:none]
        println("Blocks: $n_blocks")
        println("Draws per block: $n_sim")        
    end

    # For n_sim*n_times iterations within each block, generate a new parameter draw.
    # Decide to accept or reject, and save every (n_times)th draw that is accepted.

    all_rejections = 0

    # Initialize matrices for parameter draws and transition matrices
    para_sim = zeros(n_sim, num_parameters(m))
    like_sim = zeros(n_sim)
    post_sim = zeros(n_sim)
    TTT_sim  = zeros(n_sim, num_states_augmented(m)^2)
    RRR_sim  = zeros(n_sim, num_states_augmented(m)*num_shocks_exogenous(m))
    CCC_sim  = zeros(n_sim, num_states_augmented(m))
    z_sim    = zeros(n_sim, num_states_augmented(m))

    # Open HDF5 file for saving output
    simfile = h5open(rawpath(m,"estimate","sim_save.h5"),"w")
    n_saved_obs = n_sim * (n_blocks - n_burn)
    parasim = d_create(simfile, "parasim", datatype(Float32),
                       dataspace(n_saved_obs,n_params), "chunk", (n_sim,n_params))
    postsim = d_create(simfile, "postsim", datatype(Float32),
                       dataspace(n_saved_obs,1), "chunk", (n_sim,1))
    TTTsim  = d_create(simfile, "TTTsim", datatype(Float32),
                       dataspace(n_saved_obs,num_states_augmented(m)^2),"chunk",(n_sim,num_states_augmented(m)^2))
    RRRsim  = d_create(simfile, "RRRsim", datatype(Float32),
                       dataspace(n_saved_obs,num_states_augmented(m)*num_shocks_exogenous(m)),"chunk",
                       (n_sim,num_states_augmented(m)*num_shocks_exogenous(m)))
    zsim    = d_create(simfile, "zsim", datatype(Float32),
                       dataspace(n_saved_obs,num_states_augmented(m)),"chunk",(n_sim,num_states_augmented(m)))

    # likesim = d_create(simfile, "likesim", datatype(Float32),
    #                  dataspace(n_saved_obs,1), "chunk", (n_sim,1))
    # CCCsim  = d_create(simfile, "CCCsim", datatype(Float32),
    #                  dataspace(n_saved_obs,num_states_augmented(m)),"chunk",(n_sim,num_states_augmented(m)))



    # keep track of how long metropolis_hastings has been sampling
    total_sampling_time = 0
    
    for i = 1:n_blocks

        tic()
        
        block_rejections = 0

        for j = 1:(n_sim*n_times)

            # Draw para_new from the proposal distribution
            para_new = rand(propdist, m; cc=cc)

            # Solve the model, check that parameters are within bounds, gensys returns a
            # meaningful system, and evaluate the posterior.

            post_out = posterior!(m, para_new, YY; mh=true)
            post_new, like_new, out = post_out[:post], post_out[:like], post_out[:mats]
            
            if VERBOSITY[verbose] >= VERBOSITY[:high] 
                println("Block $i, Iteration $j: posterior = $post_new")
            end

            # Choose to accept or reject the new parameter by calculating the
            # ratio (r) of the new posterior value relative to the old one
            # We compare min(1, r) to a number drawn randomly from a
            # uniform (0, 1) distribution. This allows us to always accept
            # the new draw if its posterior value is greater than the previous draw's,
            # but it gives some probability to accepting a draw with a smaller posterior value,
            # so that we may explore tails and other local modes.

            posterior_ratio = exp(post_new - post_old)

            x = rand(m.rng)
                        
            if x < min(1.0, posterior_ratio)
                # Accept proposed jump
                para_old = para_new
                post_old = post_new
                like_old = like_new
                propdist.μ = para_new

                TTT_old = out[:TTT]
                RRR_old = out[:RRR]
                CCC_old = out[:CCC]

                zend_old = out[:zend]
                ZZ_old = out[:ZZ]
                DD_old = out[:DD]
                QQ_old = out[:QQ]

                if VERBOSITY[verbose] >= VERBOSITY[:high] 
                    println("Block $i, Iteration $j: accept proposed jump")
                end

            else
                # Reject proposed jump
                block_rejections += 1
                
                if VERBOSITY[verbose] >= VERBOSITY[:high] 
                    println("Block $i, Iteration $j: reject proposed jump")
                end
                
            end

            # Save every (n_times)th draw

            if j % n_times == 0
                draw_index = round(Int,j/n_times)
                
                like_sim[draw_index] = like_old
                post_sim[draw_index] = post_old
                para_sim[draw_index, :] = para_old'
                TTT_sim[draw_index, :] = vec(TTT_old)'
                RRR_sim[draw_index, :] = vec(RRR_old)'
                CCC_sim[draw_index, :] = vec(CCC_old)'
                z_sim[draw_index, :] = vec(zend_old)'
            end
        end

        all_rejections += block_rejections
        block_rejection_rate = block_rejections/(n_sim*n_times)

        

        ## Once every iblock times, write parameters to a file

        # Calculate starting and ending indices for this block (corresponds to a new chunk in memory)
        block_start = n_sim*(i-n_burn-1)+1
        block_end   = block_start+n_sim-1


        # Write data to file if we're past n_burn blocks
        if i > n_burn
            parasim[block_start:block_end, :]   = @compat(map(Float32,para_sim))
            postsim[block_start:block_end, :]   = @compat(map(Float32, post_sim))
            # likesim[block_start:block_end, :] = @compat(map(Float32, like_sim))
            TTTsim[block_start:block_end,:]     = @compat(map(Float32,TTT_sim))
            RRRsim[block_start:block_end,:]     = @compat(map(Float32, RRR_sim))
            zsim[block_start:block_end,:]       = @compat(map(Float32, z_sim))
        end


        block_time = toq()

        # Print status
        if VERBOSITY[verbose] > VERBOSITY[:none]

            # Calculate time to complete this block, average block
            # time, and expected time to completion

            total_sampling_time += block_time
            expected_time_remaining_sec = (total_sampling_time/i)*(n_blocks - i)
            expected_time_remaining_hrs = expected_time_remaining_sec/3600

            println("Completed $i of $n_blocks blocks.")
            println("Total time to compute $i blocks: $total_sampling_time")
            println("Expected time remaining for Metropolis-Hastings: $expected_time_remaining_hrs hours")
            println("Block $i rejection rate: $block_rejection_rate \n")
        end
        
    end # of block

    close(simfile)

    rejection_rate = all_rejections/(n_blocks*n_sim*n_times)
    if VERBOSITY[verbose] > VERBOSITY[:none]
        println("Overall rejection rate: $rejection_rate")
    end
end # of loop over blocks


#=
doc"""
compute_parameter_covariance{T<:AbstractDSGEModel}(m::T)

### Parameters
* `m::AbstractDSGEModel`: the model object

### Description:
Calculates the parameter covariance matrix from saved parameter draws, and writes it to the
parameter_covariance.h5 file in the `workpath(m)` directory.
"""
=#
function compute_parameter_covariance{T<:AbstractDSGEModel}(m::T)

    # Read in saved parameter draws
    param_draws_path = rawpath(m,"estimate","sim_save.h5")
    if !isfile(param_draws_path)
        @printf STDERR "Saved parameter draws not found."
        return
    end
    sim_h5 = h5open(param_draws_path, "r")
    param_draws = read(sim_h5, "parasim")
    close(sim_h5)
    
    # Calculate covariance matrix
    param_covariance = cov(param_draws)

    # Write to file
    cov_h5 = h5open(workpath(m, "estimate","parameter_covariance.h5"),"w")
    cov_h5["param_covariance"] = param_covariance
    close(cov_h5)

    return param_covariance;
end
