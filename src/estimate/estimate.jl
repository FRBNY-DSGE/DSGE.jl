"""
```
estimate(m, data; verbose = :low, proposal_covariance = Matrix(undef, 0, 0))
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
- `mle`: Set to true if parameters should be estimated by maximum likelihood directly.
    If this is set to true, this function will return after estimating parameters.
- `run_MH`: Set to false to disable sampling from the posterior.
"""
function estimate(m::AbstractModel, df::DataFrame;
                  verbose::Symbol = :low,
                  proposal_covariance::Matrix = Matrix(undef, 0, 0),
                  mle::Bool = false,
                  run_MH::Bool = true)
    data = df_to_matrix(m, df)
    estimate(m, data; verbose = verbose, proposal_covariance = proposal_covariance,
             mle = mle, run_MH = run_MH)
end
function estimate(m::AbstractModel;
                  verbose::Symbol = :low,
                  proposal_covariance::Matrix = Matrix(undef, 0, 0),
                  mle::Bool = false,
                  run_MH::Bool = true)
    # Load data
    df = load_data(m; verbose = verbose)
    estimate(m, df; verbose = verbose, proposal_covariance = proposal_covariance,
             mle = mle, run_MH = run_MH)
end
function estimate(m::AbstractModel, data::AbstractArray;
                  verbose::Symbol = :low,
                  proposal_covariance::Matrix = Matrix(undef, 0,0),
                  mle::Bool = false,
                  run_MH::Bool = true)

    ########################################################################################
    ### Step 1: Find posterior/likelihood mode (if reoptimizing, run optimization routine)
    ########################################################################################

    # Specify starting mode

    vint = get_setting(m, :data_vintage)
    if reoptimize(m)
        println("Reoptimizing...")

        # Inputs to optimization algorithm
        n_iterations       = get_setting(m, :optimization_iterations)
        ftol               = get_setting(m, :optimization_ftol)
        xtol               = get_setting(m, :optimization_xtol)
        gtol               = get_setting(m, :optimization_gtol)
        step_size          = get_setting(m, :optimization_step_size)
        converged          = false

        # If the algorithm stops only because we have exceeded the maximum number of
        # iterations, continue improving guess of modal parameters
        total_iterations = 0
        optimization_time = 0
        max_attempts = get_setting(m, :optimization_attempts)
        attempts = 1

        while !converged
            begin_time = time_ns()
            out, H = optimize!(m, data;
                               method = get_setting(m, :optimization_method),
                               ftol = ftol, grtol = gtol, xtol = xtol,
                               iterations = n_iterations, show_trace = true, step_size = step_size,
                               verbose = verbose,
                               mle = mle)

            attempts += 1
            total_iterations += out.iterations
            converged = out.converged || attempts > max_attempts

            end_time = (time_ns() - begin_time)/1e9

            if VERBOSITY[verbose] >= VERBOSITY[:low]
                @printf "Total iterations completed: %d\n" total_iterations
                @printf "Optimization time elapsed: %5.2f\n" optimization_time += end_time
            end

            # Write params to file after every `n_iterations` iterations
            params = map(θ->θ.value, m.parameters)
            h5open(rawpath(m, "estimate", "paramsmode.h5"),"w") do file
                file["params"] = params
            end
        end

        # write parameters to file one last time so we have the final mode
        h5open(rawpath(m, "estimate", "paramsmode.h5"),"w") do file
            file["params"] = params
        end
    end

    params = map(θ->θ.value, m.parameters)

    # Sampling does not make sense if mle=true
    if mle || !run_MH
        return nothing
    end

    ########################################################################################
    ### Step 2: Compute proposal distribution
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
	end

	# Compute inverse hessian and create proposal distribution, or
        # just create it with the given cov matrix if we have it
        propdist = if isempty(proposal_covariance)
            # Make sure the mode and hessian have the same number of parameters
            n = length(params)
            @assert (n, n) == size(hessian)

            # Compute the inverse of the Hessian via eigenvalue decomposition
            S_diag, U = eigen(hessian)
            big_eig_vals = findall(x -> x > 1e-6, S_diag)
            hessian_rank = length(big_eig_vals)

            S_inv = zeros(n, n)
            for i = (n-hessian_rank+1):n
                S_inv[i, i] = 1/S_diag[i]
            end

            hessian_inv = U*sqrt.(S_inv) #this is the inverse of the hessian
            DegenerateMvNormal(params, hessian_inv)
        else
            DegenerateMvNormal(params, proposal_covariance)
        end

        if rank(propdist) != n_parameters_free(m)
            println("problem –    shutting down dimensions")
        end

    ########################################################################################
    ### Step 3: Sample from posterior using Metropolis-Hastings algorithm
    ########################################################################################

    # Set the jump size for sampling
    cc0 = get_setting(m, :mh_cc0)
    cc = get_setting(m, :mh_cc)

    metropolis_hastings(propdist, m, data, cc0, cc; verbose=verbose)

    ########################################################################################
    ### Step 4: Calculate and save parameter covariance matrix
    ########################################################################################

    compute_parameter_covariance(m)

    return nothing
end


"""
```
metropolis_hastings(propdist::Distribution, m::AbstractModel,
    data::Matrix{T}, cc0::T, cc::T; verbose::Symbol = :low) where {T<:AbstractFloat}
```

Implements the Metropolis-Hastings MCMC algorithm for sampling from the posterior
distribution of the parameters.

### Arguments

- `propdist`: The proposal distribution that Metropolis-Hastings begins sampling from.
- `m`: The model object
- `data`: Data matrix for observables
- `cc0`: Jump size for initializing Metropolis-Hastings.
- `cc`: Jump size for the rest of Metropolis-Hastings.

### Optional Arguments

- `verbose`: The desired frequency of function progress messages printed to
  standard out. One of:

```
   - `:none`: No status updates will be reported.
   - `:low`: Status updates provided at each block.
   - `:high`: Status updates provided at each draw.
```
"""
function metropolis_hastings(propdist::Distribution,
                             m::AbstractModel,
                             data::AbstractArray,
                             cc0::T,
                             cc::T;
                             verbose::Symbol=:low) where {T<:AbstractFloat}


    # If testing, set the random seeds at fixed numbers
    if m.testing
        Random.seed!(m.rng, 654)
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

    n_blocks = n_mh_blocks(m)
    n_sim    = n_mh_simulations(m)
    n_burn   = n_mh_burn(m)
    mhthin   = mh_thin(m)

    initialized = false
    while !initialized
        post_old = posterior!(m, para_old, data; mh = true)
        if post_old > -Inf
            propdist.μ = para_old
            initialized = true
        else
            para_old = rand(propdist, m; cc=cc0)
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

    # Open HDF5 file for saving parameter draws
    simfile = h5open(rawpath(m,"estimate","mhsave.h5"),"w")
    n_saved_obs = n_sim * (n_blocks - n_burn)
    parasim = d_create(simfile, "mhparams", datatype(Float32),
                       dataspace(n_saved_obs,n_params),
                       "chunk", (n_sim,n_params))

    # Keep track of how long metropolis_hastings has been sampling
    total_sampling_time = 0.

    for block = 1:n_blocks

        begin_time = time_ns()
        block_rejections = 0

        for j = 1:(n_sim*mhthin)

            # Draw para_new from the proposal distribution
            para_new = rand(propdist, m; cc = cc)

            # Solve the model (checking that parameters are within bounds and
            # gensys returns a meaningful system) and evaluate the posterior
            post_new = posterior!(m, para_new, data; mh = true)

            if VERBOSITY[verbose] >= VERBOSITY[:high]
                println("Block $block, Iteration $j: posterior = $post_new")
            end

            # Choose to accept or reject the new parameter by calculating the
            # ratio (r) of the new posterior value relative to the old one We
            # compare min(1, r) to a number drawn randomly from a uniform (0, 1)
            # distribution. This allows us to always accept the new draw if its
            # posterior value is greater than the previous draw's, but it gives
            # some probability to accepting a draw with a smaller posterior
            # value, so that we may explore tails and other local modes.
            r = exp(post_new - post_old)
            x = rand(m.rng)

            if x < min(1.0, r)
                # Accept proposed jump
                para_old = para_new
                post_old = post_new
                propdist.μ = para_new

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
                draw_index = convert(Int, j/mhthin)
                mhparams[draw_index, :]  = para_old'
            end
        end # of block

        all_rejections += block_rejections
        block_rejection_rate = block_rejections/(n_sim*mhthin)

        ## Once every iblock times, write parameters to a file

        # Calculate starting and ending indices for this block (corresponds to a new chunk in memory)
        block_start = n_sim*(block-n_burn-1)+1
        block_end   = block_start+n_sim-1

        # Write parameters to file if we're past n_burn blocks
        if block > n_burn
            parasim[block_start:block_end, :] = map(Float32, mhparams)
        end

        # Calculate time to complete this block, average block time, and
        # expected time to completion
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            end_time = (time_ns() - begin_time)/1e9
            total_sampling_time += end_time
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

    rejection_rate = all_rejections / (n_blocks*n_sim*mhthin)
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
        @printf stderror "Saved parameter draws not found.\n"
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

"""
```
get_estimation_output_files(m)
```

Returns a `Dict{Symbol, String}` with all files created by `estimate(m)`.
"""
function get_estimation_output_files(m::AbstractModel)
    output_files = Dict{Symbol, String}()

    for file in [:paramsmode, :hessian, :mhsave]
        output_files[file] = rawpath(m, "estimate", "$file.h5")
    end

    for file in [:paramsmean, :parameter_covariance]
        output_files[file] = workpath(m, "estimate", "$file.h5")
    end

    for file in [:priors, :prior_posterior_means, :moments]
        output_files[file] = tablespath(m, "estimate", "$file.tex")
    end

    return output_files
end
