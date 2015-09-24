# This program produces and saves draws from the posterior distribution of
# the parameters.

using HDF5, Compat
include("../../test/util.jl")
doc"""
estimate{T<:AbstractDSGEModel}(m::T; verbose::Symbol=:low, testing::Bool=false)

### Parameters:
- `m`: The model object

### Optional Arguments:
- `verbose`: The desired frequency of function progress messages printed to standard out.

   - `:none`: No status updates will be reported.

   - `:low`: Status updates will be provided in csminwel and at each block in Metropolis-Hastings.

   - `:high`: Status updates provided at each iteration in Metropolis-Hastings.

- `testing`: Run `estimate()` in testing mode. In this case, a set of predetermined random numbers are read in and used rather than drawn from the Random Number Generator, and Metropolis-Hastings runs for `num_mh_simulations_test`*`num_mh_blocks_test` simulations.
"""
function estimate{T<:AbstractDSGEModel}(m::T; verbose::Symbol=:low, testing::Bool=false, using_matlab_sigscale::Bool=false)

    ###################################################################################################
    ### Step 1: Initialize
    ###################################################################################################

    # Set up levels of verbose-ness
    verboseness = Dict{Symbol,Int}(:none => 0, :low => 1, :high => 2)
    
    # Load data
    in_path = inpath(m)
    out_path = outpath(m)

    h5 = h5open(joinpath(in_path,"YY.h5"), "r") 
    YY = read(h5["YY"])
    close(h5)

    post = posterior(m, YY)


    ###################################################################################################
    ### Step 2: Find posterior mode (if reoptimizing, run csminwel)
    ###################################################################################################
    
    # Specify starting mode

    if verboseness[verbose] > verboseness[:none] 
        println("Reading in previous mode")
    end
    
    mode = []

    if m.reoptimize
        h5 = h5open("$in_path/mode_in.h5","r") 
        mode = read(h5["params"])   #it's mode in mode_in_optimized, but params in mode_in
        close(h5)
    else
        h5 = h5open("$in_path/mode_in_optimized.h5","r") 
        mode = read(h5["mode"])
        close(h5)
    end


    update!(m, mode)

    if m.reoptimize
        println("Reoptimizing...")
        
        # Inputs to minimization algorithm
        function posterior_min!{T<:AbstractFloat}(x::Vector{T})
            tomodel!(m,x)
            return -posterior(m, YY, catchGensysErrors=true)
        end

        xh = toreal(m.parameters)
        H = 1e-4 * eye(num_parameters(m))
        nit = 1000
        crit = 1e-10
        converged = false

        # If the algorithm stops only because we have exceeded the maximum number of
        # iterations, continue improving guess of modal parameters
        while !converged
            verbose_bool = verboseness[verbose] > verboseness[:none] ? true : false
            out, H = csminwel(posterior_min!, xh, H; ftol=crit, iterations=nit, show_trace=true, verbose=verbose_bool)
            xh = out.minimum
            converged = !out.iteration_converged
        end

        # Transform modal parameters so they are no longer bounded (i.e., allowed
        # to lie anywhere on the real line).
        tomodel!(m, xh)
        mode = [α.value for α in m.parameters]

        # Write mode to file
        h5open("$out_path/mode_out.h5","w") do h5
            h5["mode"] = mode
        end
        
    end

    
    ###################################################################################################
    ### Step 3: Compute proposal distribution
    ###################################################################################################

    hessian = []
    
    # Calculate the Hessian at the posterior mode
    if m.recalculate_hessian
        if verboseness[verbose] > verboseness[:none] 
            println("Recalculating Hessian...")
        end
        
        hessian, _ = hessizero!(m, mode, YY; verbose=true)

        h5open("$out_path/hessian.h5","w") do h5
            h5["hessian"] = hessian
        end

    else
        if verboseness[verbose] > verboseness[:none]
            println("Using pre-calculated Hessian")
        end

        h5 = h5open("$in_path/hessian_optimized.h5","r") 
        hessian = read(h5["hessian"])
        close(h5)

    end

    # The hessian is used to calculate the variance of the proposal
    # distribution, which is used to draw a new parameter in each iteration of
    # the algorithm.

    propdist = proposal_distribution(mode, hessian,
                                     use_matlab_sigscale=using_matlab_sigscale,
                                     sigscalepath=inpath(m), verbose=verbose)

    if propdist.rank != num_parameters_free(m)
        println("problem –    shutting down dimensions")
    end


    ###################################################################################################
    ### Step 4: Sample from posterior using Metropolis-Hastings algorithm
    ###################################################################################################
    
    # Set the jump size for sampling
    cc0 = 0.01
    cc = 0.09

    if !testing
        metropolis_hastings(propdist, m, YY, cc0, cc; verbose=verbose)
    else
        randvecs = []
        randvals = []
        
        h5= h5open("/data/dsge_data_dir/dsgejl/estimate/save/input_data/rand_save_big.h5","r")
        randvecs = read(h5["randvecs"])
        randvals = read(h5["randvals"])
        close(h5)

        metropolis_hastings(propdist, m, YY, cc0, cc; verbose=verbose, randvecs=randvecs, randvals=randvals,
                            use_matlab_sigscale=using_matlab_sigscale)
        
    end

    
    ###################################################################################################
    ### Step 5: Calculate and save parameter covariance matrix
    ###################################################################################################
    
    # Set up HDF5 file for saving
    h5path = joinpath(out_path,"sim_save.h5")

    # Read in saved parameter draws
    sim_h5 = h5open(h5path, "r+")

    θ = read(sim_h5, "parasim")

    # Calculate covariance matrix
    cov_θ = cov(θ)
    write(sim_h5, "cov_θ", map(Float32,cov_θ))   #Save as single-precision float matrix

    # Close the file
    close(sim_h5)
end

doc"""
proposal_distribution{T<:AbstractFloat, V<:AbstractString}(μ::Vector{T}, hessian::Matrix{T}; use_matlab_sigscale::Bool=false, sigscalepath::V="", verbose=:low)

### Parameters:
* `μ`: Vector of means
* `h` Hessian matrix

### Description:
Compute proposal distribution: degenerate normal with mean μ and covariance hessian^(-1)
"""
function proposal_distribution{T<:AbstractFloat, V<:AbstractString}(μ::Vector{T}, hessian::Matrix{T}; use_matlab_sigscale::Bool=false, sigscalepath::V="", verbose=:low)

    # Set up levels of verbose-ness
    verboseness = Dict{Symbol,Int}(:none => 0, :low => 1, :high => 2)
    
    n = length(μ)
    @assert (n, n) == size(hessian)

    S_diag, U = eig(hessian)
    big_evals = find(x -> x > 1e-6, S_diag)
    rank = length(big_evals)
    
    S_inv = zeros(n, n)
    for i = (n-rank+1):n
        S_inv[i, i] = 1/S_diag[i]
    end

    σ = U*sqrt(S_inv)

    # Use predefined reference sigscale if indicated
    if use_matlab_sigscale

        if verboseness[verbose] > verboseness[:none]
            println("Using sigscale from Matlab")
        end

        h5 = h5open(joinpath(sigscalepath,"sigscale.h5"), "r")
        σ = read(h5, "sigscale")
        close(h5)

    end

    return DegenerateMvNormal(μ, σ, rank)
end

doc"""
metropolis_hastings{T<:AbstractFloat}(propdist::Distribution, m::AbstractDSGEModel, YY::Matrix{T}, cc0::T, cc::T; randvecs = [], randvals = [], verbose = :low)

### Parameters
* `propdist` The proposal distribution that Metropolis-Hastings begins sampling from.
* `m`: The model object
* `YY`: Data matrix for observables
* `cc0`: Jump size for initializing Metropolis-Hastings.
* `cc`: Jump size for the rest of Metropolis-Hastings.

### Optional Arguments
* `randvecs`: A set of predetermined random vectors drawn from N(0,I) used to avoid true randomness during testing. 
* `randvals`: A set of predetermined random values drawn from U(0,1) used during testing.
* `verbose`: The desired frequency of function progress messages printed to standard out.

   - `:none`: No status updates will be reported.

   - `:low`: Status updates provided at each block.

   - `:high`: Status updates provided at each draw.
 
### Description
Implements the Metropolis-Hastings MCMC algorithm.
"""
function metropolis_hastings{T<:AbstractFloat}(propdist::Distribution, m::AbstractDSGEModel,
    YY::Matrix{T}, cc0::T, cc::T; randvecs = [], randvals = [], verbose = :low, use_matlab_sigscale=false)

    # Set up levels of verbose-ness
    verboseness = Dict{Symbol,Int}(:none => 0, :low => 1, :high => 2)
    
    # If testing, then we read in a specific sequence of "random" vectors and numbers
    testing = !(randvecs == [] && randvals == [])

    if verboseness[verbose] > verboseness[:none]
        println("Testing = $testing")
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
    para_old = propdist.μ
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
        if testing
            para_old = propdist.μ + cc0*propdist.σ*randvecs[:, 1]
            n_blocks = m.num_mh_blocks_test
            n_sim = m.num_mh_simulations_test
            n_burn = m.num_mh_burn_test
            n_times = m.mh_thinning_step
        else
            para_old = rand(propdist; cc=cc0)
            n_blocks = m.num_mh_blocks
            n_sim = m.num_mh_simulations
            n_burn = m.num_mh_burn
            n_times = m.mh_thinning_step
        end

        post_old, like_old, out = posterior!(m, para_old, YY; mh=true)
        
        if post_old > -Inf
            propdist.μ = para_old

            TTT_old = out[:TTT]
            RRR_old = out[:RRR]
            CCC_old = out[:CCC]
            zend_old  = out[:zend]

            initialized = true
        end

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

    # # Open HDF5 file for saving output
    
    h5path = joinpath(outpath(m),"sim_save.h5")

    simfile = h5open(h5path,"w")

    n_saved_obs = n_sim * (n_blocks - n_burn)

    parasim = d_create(simfile, "parasim", datatype(Float32),
                       dataspace(n_saved_obs,n_params), "chunk", (n_sim,n_params))

    # likesim = d_create(simfile, "likesim", datatype(Float32),
    #                  dataspace(n_saved_obs,1), "chunk", (n_sim,1))

    postsim = d_create(simfile, "postsim", datatype(Float32),
                       dataspace(n_saved_obs,1), "chunk", (n_sim,1))

    TTTsim  = d_create(simfile, "TTTsim", datatype(Float32),
                       dataspace(n_saved_obs,num_states_augmented(m)^2),"chunk",(n_sim,num_states_augmented(m)^2))

    RRRsim  = d_create(simfile, "RRRsim", datatype(Float32),
                       dataspace(n_saved_obs,num_states_augmented(m)*num_shocks_exogenous(m)),"chunk",
                       (n_sim,num_states_augmented(m)*num_shocks_exogenous(m)))

    # CCCsim  = d_create(simfile, "CCCsim", datatype(Float32),
    #                  dataspace(n_saved_obs,num_states_augmented(m)),"chunk",(n_sim,num_states_augmented(m)))

    zsim    = d_create(simfile, "zsim", datatype(Float32),
                       dataspace(n_saved_obs,num_states_augmented(m)),"chunk",(n_sim,num_states_augmented(m)))

    if testing
        rows, cols = size(randvecs)
        numvals = size(randvals)[1]
    end

    # keep track of how long metropolis_hastings has been sampling
    total_sampling_time = 0
    
    for i = 1:n_blocks

        tic()
        
        block_rejections = 0

        for j = 1:(n_sim*n_times)

            # Draw para_new from the proposal distribution

            if testing
                para_new = propdist.μ + cc*propdist.σ*randvecs[:, mod(j,cols)+1]
            else
                para_new = rand(propdist; cc=cc)
            end

            # Solve the model, check that parameters are within bounds, gensys returns a
            # meaningful system, and evaluate the posterior.

            post_new, like_new, out = posterior!(m, para_new, YY; mh=true)
            
            if verboseness[verbose] >= verboseness[:high] 
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

            if testing
                k = (i-1)*(n_sim*n_times) + mod(j,numvals)+1
                x = randvals[mod(j,numvals)+1]
            else
                x = rand(m.rng)
            end

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

                if verboseness[verbose] >= verboseness[:high] 
                    println("Block $i, Iteration $j: accept proposed jump")
                end

            else
                # Reject proposed jump
                block_rejections += 1
                
                if verboseness[verbose] >= verboseness[:high] 
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
            parasim[block_start:block_end, :] = map(Float32,para_sim)
            postsim[block_start:block_end, :] = map(Float32, post_sim)
            # likesim[block_start:block_end, :] = map(Float32, like_sim)
            TTTsim[block_start:block_end,:]  = map(Float32,TTT_sim)
            RRRsim[block_start:block_end,:]  = map(Float32, RRR_sim)
            zsim[block_start:block_end,:]  = map(Float32, z_sim)
        end


        block_time = toq()

        # Print status
        if verboseness[verbose] > verboseness[:none]

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
    if verboseness[verbose] > verboseness[:none]
        println("Overall rejection rate: $rejection_rate")
    end
end # of loop over blocks
