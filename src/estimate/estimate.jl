# This program produces and saves draws from the posterior distribution of
# the parameters.

"""
```
estimate(m::AbstractModel; verbose::Symbol=:low, proposal_covariance=Matrix())
```

This routine implements the full estimation stage of the FRBNY DSGE model.

### Arguments
- `m`: The model object

### Optional Arguments:
- `verbose`: The desired frequency of function progress messages printed to standard out.

   - `:none`: No status updates will be reported.

   - `:low`: Status updates will be provided in csminwel and at each block in Metropolis-Hastings.

   - `:high`: Status updates provided at each iteration in Metropolis-Hastings.

- `proposal_covariance`: Used to test the metropolis_hastings algorithm with a precomputed
  covariance matrix for the proposal distribution. When the Hessian is singular, eigenvectors
  corresponding to zero eigenvectors are not well defined, so eigenvalue decomposition can
  cause problems. Passing a precomputed matrix allows us to ensure that the rest of the
  routine has not broken.
"""
function estimate(m::AbstractModel; 
                  verbose::Symbol=:low,
                  proposal_covariance::Matrix=Matrix())

    ########################################################################################
    ### Step 1: Initialize
    ########################################################################################

    # Load data
    YY = h5open(inpath(m, "data","data_$(data_vintage(m)).h5"), "r") do f
        read(f, "YY")
    end

    post = posterior(m, YY)[:post]


    ########################################################################################
    ### Step 2: Find posterior mode (if reoptimizing, run csminwel)
    ########################################################################################
    
    # Specify starting mode

    if VERBOSITY[verbose] >= VERBOSITY[:low] 
        println("Reading in previous mode")
    end
    
    params = if optimize(m)
        #it's mode in params_mode, but params in params_start
        h5open(inpath(m, "user", "params_start.h5"),"r") do file
            read(file, "params")   
        end

    else
        h5open(inpath(m, "user", "params_mode.h5"),"r") do file
            read(file, "mode")
        end
    end

    update!(m, params)

    if optimize(m)
        println("Reoptimizing...")
        
        # Inputs to optimization algorithm
        n_iterations       = 100
        ftol      = 1e-10
        converged = false

        # If the algorithm stops only because we have exceeded the maximum number of
        # iterations, continue improving guess of modal parameters
        total_iterations = 0
        tic()
        while !converged
            out, H = optimize!(m, YY; 
                ftol=ftol, iterations=n_iterations, show_trace=true, verbose=verbose)
            converged = !out.iteration_converged

            total_iterations += out.iterations
            if VERBOSITY[verbose] >= VERBOSITY[:low] 
                @printf "Total iterations completed: %d" total_iterations
                @printf "Optimization time elapsed: %5.2f" toq()
            end

            # Write params to file after every `n_iterations` iterations
            params = map(θ->θ.value, m.parameters)
            h5open(rawpath(m, "estimate", "params_mode_out.h5"),"w") do file
                file["mode"] = params
            end
        end
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
    hessian = if calculate_hessian(m)
        if VERBOSITY[verbose] >= VERBOSITY[:low] 
            println("Recalculating Hessian...")
        end
        
        hessian, _ = hessian!(m, params, YY; verbose=verbose)

        h5open(rawpath(m, "estimate","hessian.h5"),"w") do file
            file["hessian"] = hessian
        end

        hessian

    # Read in a pre-optimized mode
    else
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            println("Using pre-calculated Hessian")
        end

        hessian = h5open(inpath(m, "user", "hessian.h5"),"r") do file
            read(file, "hessian")
        end
    end

    # Compute inverse hessian and create proposal distribution, or
    # just create it with the given cov matrix if we have it
    propdist = if isempty(proposal_covariance)
        # Make sure the mode and hessian have the same number of parameters
        n = length(params)
        @assert (n, n) == size(hessian)

        # Compute the inverse of the Hessian via eigenvalue decomposition
        S_diag, U = eig(hessian)
        big_eig_vals = find(x -> x > 1e-6, S_diag)
        rank = length(big_eig_vals)
    
        S_inv = zeros(n, n)

        for i = (n-rank+1):n
            S_inv[i, i] = 1/S_diag[i]
        end
        
        hessian_inv = U*sqrt(S_inv) #this is the inverse of the hessian
        DegenerateMvNormal(params, hessian_inv, rank)
    else
        DegenerateMvNormal(params, proposal_covariance)
    end
    
    if rank(propdist) != n_parameters_free(m)
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


"""
```
metropolis_hastings{T<:AbstractFloat}(propdist::Distribution, m::AbstractModel,
    YY::Matrix{T}, cc0::T, cc::T; verbose::Symbol = :low)
```

Implements the Metropolis-Hastings MCMC algorithm for sampling from the posterior
distribution of the parameters.

### Arguments
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
"""
function metropolis_hastings{T<:AbstractFloat}(propdist::Distribution,
                                               m::AbstractModel,
                                               YY::Matrix{T},
                                               cc0::T,
                                               cc::T;
                                               verbose::Symbol=:low)

    # If testing, set the random seeds at fixed numbers
    if m.testing
        srand(m.rng, 654)
    end
        
    # Set number of draws, how many we will save, and how many we will burn
    # (initialized here for scoping; will re-initialize in the while loop)
    n_blocks = 0
    n_sim = 0
    mhthin = 0
    n_burn = 0
    n_params = n_parameters(m)

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

        n_blocks = n_mh_blocks(m)
        n_sim    = n_mh_simulations(m)
        n_burn   = n_mh_burn(m)
        mhthin  = mh_thin(m)

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
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("Blocks: $n_blocks")
        println("Draws per block: $n_sim")        
    end

    # For n_sim*mhthin iterations within each block, generate a new parameter draw.
    # Decide to accept or reject, and save every (mhthin)th draw that is accepted.

    all_rejections = 0

    # Initialize matrices for parameter draws and transition matrices
    mhparams = zeros(n_sim, n_parameters(m))
    mhlikelihood = zeros(n_sim)
    mhposterior = zeros(n_sim)
    mhTTT  = zeros(n_sim, n_states_augmented(m)^2)
    mhRRR  = zeros(n_sim, n_states_augmented(m)*n_shocks_exogenous(m))
    CCC_sim  = zeros(n_sim, n_states_augmented(m))
    mhzend    = zeros(n_sim, n_states_augmented(m))

    println("mode it to here")
    # Open HDF5 file for saving output
    simfile = h5open(rawpath(m,"estimate","mh_save.h5"),"w")
    n_saved_obs = n_sim * (n_blocks - n_burn)
    parasim = d_create(simfile, "parasim", datatype(Float32),
                       dataspace(n_saved_obs,n_params), "chunk", (n_sim,n_params))
    postsim = d_create(simfile, "postsim", datatype(Float32),
                       dataspace(n_saved_obs,1), "chunk", (n_sim,1))
    TTTsim  = d_create(simfile, "TTTsim", datatype(Float32),
                       dataspace(n_saved_obs,n_states_augmented(m)^2),"chunk",(n_sim,n_states_augmented(m)^2))
    RRRsim  = d_create(simfile, "RRRsim", datatype(Float32),
                       dataspace(n_saved_obs,n_states_augmented(m)*n_shocks_exogenous(m)),"chunk",
                       (n_sim,n_states_augmented(m)*n_shocks_exogenous(m)))
    zsim    = d_create(simfile, "zsim", datatype(Float32),
                       dataspace(n_saved_obs,n_states_augmented(m)),"chunk",(n_sim,n_states_augmented(m)))

    # likesim = d_create(simfile, "likesim", datatype(Float32),
    #                  dataspace(n_saved_obs,1), "chunk", (n_sim,1))
    # CCCsim  = d_create(simfile, "CCCsim", datatype(Float32),
    #                  dataspace(n_saved_obs,n_states_augmented(m)),"chunk",(n_sim,n_states_augmented(m)))

    # keep track of how long metropolis_hastings has been sampling
    total_sampling_time = 0
    
    for block = 1:n_blocks

        tic()
        block_rejections = 0

        for j = 1:(n_sim*mhthin)

            # Draw para_new from the proposal distribution
            para_new = rand(propdist, m; cc=cc)

            # Solves the model, check that parameters are within bounds, gensys returns a
            # meaningful system, and evaluate the posterior.
            post_out = posterior!(m, para_new, YY; mh=true)
            post_new, like_new, out = post_out[:post], post_out[:like], post_out[:mats]
            
            if VERBOSITY[verbose] >= VERBOSITY[:high] 
                println("Block $block, Iteration $j: posterior = $post_new")
            end

            # Choose to accept or reject the new parameter by calculating the
            # ratio (r) of the new posterior value relative to the old one
            # We compare min(1, r) to a number drawn randomly from a
            # uniform (0, 1) distribution. This allows us to always accept
            # the new draw if its posterior value is greater than the previous draw's,
            # but it gives some probability to accepting a draw with a smaller posterior value,
            # so that we may explore tails and other local modes.

            r = exp(post_new - post_old)

            x = rand(m.rng)
                        
            if x < min(1.0, r)
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
                    println("Block $block, Iteration $j: accept proposed jump")
                end

            else
                # Reject proposed jump
                block_rejections += 1
                
                if VERBOSITY[verbose] >= VERBOSITY[:high] 
                    println("Block $block, Iteration $j: reject proposed jump")
                end
                
            end

            # Save every (mhthin)th draw

            if j % mhthin == 0
                draw_index = round(Int,j/mhthin)
                
                mhlikelihood[draw_index] = like_old
                mhposterior[draw_index] = post_old
                mhparams[draw_index, :] = para_old'
                mhTTT[draw_index, :] = vec(TTT_old)'
                mhRRR[draw_index, :] = vec(RRR_old)'
                CCC_sim[draw_index, :] = vec(CCC_old)'
                mhzend[draw_index, :] = vec(zend_old)'
            end
        end

        all_rejections += block_rejections
        block_rejection_rate = block_rejections/(n_sim*mhthin)

        ## Once every iblock times, write parameters to a file

        # Calculate starting and ending indices for this block (corresponds to a new chunk in memory)
        block_start = n_sim*(block-n_burn-1)+1
        block_end   = block_start+n_sim-1

        # Write data to file if we're past n_burn blocks
        if i > n_burn
            parasim[block_start:block_end, :]   = @compat(map(Float32,mhparams))
            postsim[block_start:block_end, :]   = @compat(map(Float32, mhposterior))
            # likesim[block_start:block_end, :] = @compat(map(Float32, mhlikelihood))
            TTTsim[block_start:block_end,:]     = @compat(map(Float32,mhTTT))
            RRRsim[block_start:block_end,:]     = @compat(map(Float32, mhRRR))
            zsim[block_start:block_end,:]       = @compat(map(Float32, mhzend))
        end


        # Calculate time to complete this block, average block
        # time, and expected time to completion
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            block_time = toq()
            total_sampling_time += block_time
            expected_time_remaining_sec = (total_sampling_time/block)*(n_blocks - block)
            expected_time_remaining_hrs = expected_time_remaining_sec/3600

            println("Completed $block of $n_blocks blocks.")
            println("Total time to compute $block blocks: $total_sampling_time")
            println("Expected time remaining for Metropolis-Hastings: $expected_time_remaining_hrs hours")
            println("Block $block rejection rate: $block_rejection_rate \n")
        end
        
    end # of block

    close(simfile)

    rejection_rate = all_rejections/(n_blocks*n_sim*mhthin)
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("Overall rejection rate: $rejection_rate")
    end
end # of loop over blocks

"""
```
compute_parameter_covariance{T<:AbstractModel}(m::T)
```

Calculates the parameter covariance matrix from saved parameter draws, and writes it to the
parameter_covariance.h5 file in the `workpath(m)` directory.

### Arguments
* `m::AbstractModel`: the model object
"""
function compute_parameter_covariance{T<:AbstractModel}(m::T)

    # Read in saved parameter draws
    param_draws_path = rawpath(m,"estimate","mh_save.h5")
    if !isfile(param_draws_path)
        @printf STDERR "Saved parameter draws not found."
        return
    end
    param_draws = h5open(param_draws_path, "r") do f
        read(f, "parasim")
    end
    
    # Calculate covariance matrix
    param_covariance = cov(param_draws)

    # Write to file
    h5open(workpath(m, "estimate","parameter_covariance.h5"),"w") do f
        f["param_covariance"] = param_covariance
    end

    return param_covariance;
end
