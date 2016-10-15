"""
```
estimate(m, data; verbose=:low, proposal_covariance=Matrix())
```

Estimate the DSGE parameter posterior distribution.

### Arguments:
- `m::AbstractModel`: model object

### Optional Arguments:
- `data`: well-formed data as `Matrix` or `DataFrame`. If this is not provided, the `load_data` routine will be executed.

### Keyword Arguments:
- `verbose::Symbol`: The desired frequency of function progress messages printed to standard out.
   - `:none`: No status updates will be reported.
   - `:low`: Status updates will be provided in csminwel and at each block in
     Metropolis-Hastings.
   - `:high`: Status updates provided at each iteration in Metropolis-Hastings.
- `proposal_covariance::Matrix`: Used to test the metropolis_hastings algorithm with a precomputed
  covariance matrix for the proposal distribution. When the Hessian is singular,
  eigenvectors corresponding to zero eigenvectors are not well defined, so eigenvalue
  decomposition can cause problems. Passing a precomputed matrix allows us to ensure that
  the rest of the routine has not broken. 
"""
function estimate(m::AbstractModel, df::DataFrame;
                  verbose::Symbol=:low,
                  proposal_covariance::Matrix=Matrix())
    data = df_to_matrix(m, df)
    estimate(m, data; verbose=verbose, proposal_covariance=proposal_covariance)
end
function estimate(m::AbstractModel;
                  verbose::Symbol=:low,
                  proposal_covariance::Matrix=Matrix())
    # Load data
    df = load_data(m; verbose=verbose)
    estimate(m, df; verbose=verbose, proposal_covariance=proposal_covariance)
end
function estimate(m::AbstractModel, data::Matrix{Float64};
                  verbose::Symbol=:low,
                  proposal_covariance::Matrix=Matrix())

    ########################################################################################
    ### Step 1: Initialize
    ########################################################################################

    post = posterior(m, data)[:post]

    ########################################################################################
    ### Step 2: Find posterior mode (if reoptimizing, run csminwel)
    ########################################################################################

    # Specify starting mode

    vint = get_setting(m, :data_vintage)
    if reoptimize(m)
        println("Reoptimizing...")

        # Inputs to optimization algorithm
        n_iterations       = 100
        ftol               = 1e-10
        converged          = false

        # If the algorithm stops only because we have exceeded the maximum number of
        # iterations, continue improving guess of modal parameters
        total_iterations = 0
        optimization_time = 0
        while !converged
            tic()
            out, H = optimize!(m, data;
                ftol=ftol, iterations=n_iterations, show_trace=true, verbose=verbose)
            converged = !out.iteration_converged

            total_iterations += out.iterations
            if VERBOSITY[verbose] >= VERBOSITY[:low]
                @printf "Total iterations completed: %d\n" total_iterations
                @printf "Optimization time elapsed: %5.2f\n" optimization_time += toq()
            end

            # Write params to file after every `n_iterations` iterations
            params = map(θ->θ.value, m.parameters)
            h5open(rawpath(m, "estimate", "paramsmode.h5"),"w") do file
                file["params"] = params
            end
        end
    end
    
    params = map(θ->θ.value, m.parameters)
    
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

        hessian, _ = hessian!(m, params, data; verbose=verbose)

        h5open(rawpath(m, "estimate","hessian.h5"),"w") do file
            file["hessian"] = hessian
        end

        hessian

    # Read in a pre-calculated Hessian
    else
        fn = hessian_path(m)
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            println("Using pre-calculated Hessian from $fn")
        end

        hessian = h5open(fn,"r") do file
            read(file, "hessian")
        end

        hessian
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
        DSGE.DegenerateMvNormal(params, hessian_inv)
    else
        DSGE.DegenerateMvNormal(params, proposal_covariance)
    end
    
    if DSGE.rank(propdist) != n_parameters_free(m)
        println("problem –    shutting down dimensions")
    end

    ########################################################################################
    ### Step 4: Sample from posterior using Metropolis-Hastings algorithm
    ########################################################################################

    # Set the jump size for sampling
    cc0 = 0.01
    cc = 0.09

    metropolis_hastings(propdist, m, data, cc0, cc; verbose=verbose);

    ########################################################################################
    ### Step 5: Calculate and save parameter covariance matrix
    ########################################################################################

    compute_parameter_covariance(m);

    return nothing
end


"""
```
metropolis_hastings{T<:AbstractFloat}(propdist::Distribution, m::AbstractModel,
    data::Matrix{T}, cc0::T, cc::T; verbose::Symbol = :low)
```
Implements the Metropolis-Hastings MCMC algorithm for sampling from the posterior
distribution of the parameters.

### Arguments
- `propdist` The proposal distribution that Metropolis-Hastings begins sampling from.
- `m`: The model object
- `data`: Data matrix for observables
- `cc0`: Jump size for initializing Metropolis-Hastings.
- `cc`: Jump size for the rest of Metropolis-Hastings.

### Optional Arguments
- `verbose`: The desired frequency of function progress messages printed to standard out.
   - `:none`: No status updates will be reported.
   - `:low`: Status updates provided at each block.
   - `:high`: Status updates provided at each draw.
"""
function metropolis_hastings{T<:AbstractFloat}(propdist::Distribution,
                                               m::AbstractModel,
                                               data::Matrix{T},
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
    n_sim    = 0
    mhthin   = 0
    n_burn   = 0
    n_params = n_parameters(m)

    # Initialize algorithm by drawing para_old from a normal distribution centered on the
    # posterior mode until the parameters are within bounds or the posterior value is
    # sufficiently large.
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
        mhthin   = mh_thin(m)

        post_out = posterior!(m, para_old, data; mh=true)
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
    mhparams     = zeros(n_sim, n_parameters(m))
    mhlikelihood = zeros(n_sim)
    mhposterior  = zeros(n_sim)
    mhTTT        = zeros(n_sim, n_states_augmented(m), n_states_augmented(m))
    mhRRR        = zeros(n_sim, n_states_augmented(m), n_shocks_exogenous(m))
    mhCCC        = zeros(n_sim, n_states_augmented(m))
    mhzend       = zeros(n_sim, n_states_augmented(m))

    # Open HDF5 file for saving output
    simfile = h5open(rawpath(m,"estimate","mhsave.h5"),"w")
    n_saved_obs = n_sim * (n_blocks - n_burn)
    parasim = d_create(simfile, "mhparams", datatype(Float32),
                       dataspace(n_saved_obs,n_params),
                       "chunk", (n_sim,n_params))
    postsim = d_create(simfile, "mhposterior", datatype(Float32),
                       dataspace(n_saved_obs,1),
                       "chunk", (n_sim,1))
    TTTsim  = d_create(simfile, "mhTTT", datatype(Float32),
                       dataspace(n_saved_obs,n_states_augmented(m),n_states_augmented(m)),
                       "chunk", (n_sim,n_states_augmented(m),n_states_augmented(m)))
    RRRsim  = d_create(simfile, "mhRRR", datatype(Float32),
                       dataspace(n_saved_obs,n_states_augmented(m),n_shocks_exogenous(m)),
                       "chunk", (n_sim,n_states_augmented(m),n_shocks_exogenous(m)))
    zsim    = d_create(simfile, "mhzend", datatype(Float32),
                       dataspace(n_saved_obs,n_states_augmented(m)),
                       "chunk", (n_sim,n_states_augmented(m)))

    # likesim = d_create(simfile, "mhlikelihood", datatype(Float32),
    #                  dataspace(n_saved_obs,1),
    #                  "chunk", (n_sim,1))
    # CCCsim  = d_create(simfile, "mhCCC", datatype(Float32),
    #                  dataspace(n_saved_obs,n_states_augmented(m)),
    #                  "chunk",(n_sim,n_states_augmented(m)))

    # keep track of how long metropolis_hastings has been sampling
    total_sampling_time = 0.

    for block = 1:n_blocks

        tic()
        block_rejections = 0

        for j = 1:(n_sim*mhthin)

            # Draw para_new from the proposal distribution
            para_new = rand(propdist, m; cc=cc)

            # Solves the model, check that parameters are within bounds, gensys returns a
            # meaningful system, and evaluate the posterior.
            post_out = posterior!(m, para_new, data; mh=true)
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
                draw_index = convert(Int,j/mhthin)

                mhlikelihood[draw_index] = like_old
                mhposterior[draw_index]  = post_old
                mhparams[draw_index, :]  = para_old'
                mhTTT[draw_index, :, :]  = TTT_old
                mhRRR[draw_index, :, :]  = RRR_old
                mhCCC[draw_index, :, :]  = CCC_old'
                mhzend[draw_index, :]    = zend_old'
            end
        end # of block

        all_rejections += block_rejections
        block_rejection_rate = block_rejections/(n_sim*mhthin)

        ## Once every iblock times, write parameters to a file

        # Calculate starting and ending indices for this block (corresponds to a new chunk in memory)
        block_start = n_sim*(block-n_burn-1)+1
        block_end   = block_start+n_sim-1

        # Write data to file if we're past n_burn blocks
        if block > n_burn
            parasim[block_start:block_end, :]   = map(Float32, mhparams)
            postsim[block_start:block_end, :]   = map(Float32, mhposterior)
            TTTsim[block_start:block_end,:,:]   = map(Float32, mhTTT)
            RRRsim[block_start:block_end,:,:]   = map(Float32, mhRRR)
            zsim[block_start:block_end,:]       = map(Float32, mhzend)
            # likesim[block_start:block_end, :] = map(Float32, mhlikelihood)
            # CCCsim[block_start:block_end, :]  = map(Float32, mhCCC)
        end


        # Calculate time to complete this block, average block
        # time, and expected time to completion
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            block_time = toq()
            total_sampling_time += block_time
            total_sampling_time_minutes = total_sampling_time/60
            expected_time_remaining_sec     = (total_sampling_time/block)*(n_blocks - block)
            expected_time_remaining_minutes = expected_time_remaining_sec/60

            println("Completed $block of $n_blocks blocks.")
            println("Total time to compute $block blocks: $total_sampling_time_minutes minutes")
            println("Expected time remaining for Metropolis-Hastings: $expected_time_remaining_minutes minutes")
            println("Block $block rejection rate: $block_rejection_rate \n")
        end

    end # of loop over blocks

    close(simfile)

    rejection_rate = all_rejections/(n_blocks*n_sim*mhthin)
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("Overall rejection rate: $rejection_rate")
    end
end

"""
```
compute_parameter_covariance(m::AbstractModel)
```

Calculates the parameter covariance matrix from saved parameter draws, and writes it to the
parameter_covariance.h5 file in the `workpath(m, "estimate")` directory.

### Arguments
* `m::AbstractModel`: the model object
"""
function compute_parameter_covariance(m::AbstractModel)

    # Read in saved parameter draws
    param_draws_path = rawpath(m,"estimate","mhsave.h5")
    if !isfile(param_draws_path)
        @printf STDERR "Saved parameter draws not found.\n"
        return
    end
    param_draws = h5open(param_draws_path, "r") do f
        read(f, "mhparams")
    end

    # Calculate covariance matrix
    param_covariance = cov(param_draws)

    # Write to file
    h5open(workpath(m, "estimate","parameter_covariance.h5"),"w") do f
        f["mhcov"] = param_covariance
    end

    return param_covariance
end
