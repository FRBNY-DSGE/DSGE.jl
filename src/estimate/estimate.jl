# This program produces and saves draws from the posterior distribution of
# the parameters.

using HDF5
using Debug

function estimate{T<:AbstractDSGEModel}(m::T; verbose=false)

    ### Step 1: Initialize

    # Load data
    mf = MatFile("$inpath/YY.mat")
    YY = get_variable(mf, "YY")
    close(mf)

    post = posterior(m, YY)

    ### Step 2: Find posterior mode

    # Specify starting mode
    println("Reading in previous mode...")
    mf = MatFile("$inpath/mode_in.mat")
    mode = get_variable(mf, "params")
    close(mf)
    update!(m, mode)

    if m.reoptimize
        println("Reoptimizing...")

        # Inputs to minimization algorithm
        function posterior_min!{T<:FloatingPoint}(x::Vector{T})
            tomodel!(x, m.parameters)
            return -posterior(m, YY)
        end

        xh = toreal(m.parameters)
        H = 1e-4 * eye(num_parameters(m))
        nit = 1000
        crit = 1e-10
        converged = false

        # If the algorithm stops only because we have exceeded the maximum number of
        # iterations, continue improving guess of modal parameters
        while !converged
            out, H = csminwel(posterior_min!, xh, H; ftol=crit, iterations=nit, show_trace=true, verbose=verbose)
            xh = out.minimum
            converged = !out.iteration_converged
        end

        # Transform modal parameters so they are no longer bounded (i.e., allowed
        # to lie anywhere on the real line).
        tomodel!(xh, m.parameters)
        mode = [α.value for α in m.parameters]

        # Write mode to file
        mf = MatFile("$outpath/mode_out.mat", "w")
        put_variable(mf, "mode", mode)
        close(mf)
    end

    ### Step 3: Compute proposal distribution

    # Calculate the Hessian at the posterior mode
    if m.recalculate_hessian
        println("Recalculating Hessian...")
        hessian, _ = hessizero!(m, mode, YY; verbose=true)

        mf = MatFile("$outpath/hessian.mat", "w")
        put_variable(mf, "hessian", hessian)
        close(mf)
    else
        println("Using pre-calculated Hessian...")
        mf = MatFile("$inpath/hessian.mat")
        hessian = get_variable(mf, "hessian")
        close(mf)
    end

    # The hessian is used to calculate the variance of the proposal
    # distribution, which is used to draw a new parameter in each iteration of
    # the algorithm.
    propdist = proposal_distribution(mode, hessian)
    if propdist.rank != num_parameters_free(m)
        println("problem – shutting down dimensions")
    end

    ### Step 4: Metropolis-Hastings algorithm

    # Set the jump size for sampling
    cc0 = 0.01
    cc = 0.09

    metropolis_hastings(propdist, m, YY, cc0, cc; verbose=verbose)

    # Set up HDF5 file for saving
    h5path = joinpath(outpath,"sim_save.h5")

    ### Step 5: Calculate parameter covariance matrix
    # Read in saved parameter draws
    sim_h5 = h5open(h5path, "r+")
    θ = read(sim_h5, "parasim")

    # Calculate covariance matrix
    cov_θ = cov(θ)
    write(sim_h5, "cov_θ", convert(Matrix{Float32}, cov_θ))   #Save as single-precision float matrix

    # Close the file
    close(sim_h5)
end

# Compute proposal distribution: degenerate normal with mean μ and covariance hessian^(-1)
function proposal_distribution{T<:FloatingPoint}(μ::Vector{T}, hessian::Matrix{T})
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

    return DegenerateMvNormal(μ, σ, rank)
end

function metropolis_hastings{T<:FloatingPoint}(propdist::Distribution, m::AbstractDSGEModel,
    YY::Matrix{T}, cc0::T, cc::T; randvecs = [], randvals = [], verbose = false)

    # If testing, then we read in a specific sequence of "random" vectors and numbers
    testing = !(randvecs == [] && randvals == [])
    println("Testing = $testing")

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

            TTT_old = out["TTT"]
            RRR_old = out["RRR"]
            CCC_old = out["CCC"]
            zend_old  = out["zend"]

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
    # if testing
    #     savepath  = joinpath(pwd(),"save")
    #     inpath    = savepath;
    #     outpath   = savepath;
    #     tablepath = savepath;
    #     plotpath  = savepath;
    #     logpath   = savepath;
    # end

    h5path = joinpath("$outpath","sim_save.h5")
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

    for i = 1:n_blocks
        block_rejections = 0

        for j = 1:(n_sim*n_times)

            # Draw para_new from the proposal distribution
            if testing
                para_new = propdist.μ + cc*propdist.σ*randvecs[:, mod(j,cols)]
            else
                para_new = rand(propdist; cc=cc)
            end

            # Solve the model, check that parameters are within bounds, gensys returns a
            # meaningful system, and evaluate the posterior.
            post_new, like_new, out = posterior!(m, para_new, YY; mh=true)

            if verbose 
                println("Iteration $j: posterior = $post_new")
            end

            # Choose to accept or reject the new parameter by calculating the
            # ratio (r) of the new posterior value relative to the old one
            # We compare min(1, r) to a number drawn randomly from a
            # uniform (0, 1) distribution. This allows us to always accept
            # the new draw if its posterior value is greater than the previous draw's,
            # but it gives some probability to accepting a draw with a smaller posterior value,
            # so that we may explore tails and other local modes.
            r = exp(post_new - post_old)

            if testing
                k = (i-1)*(n_sim*n_times) + j
                x = randvals[mod(j,numvals)]
            else
                x = rand()
            end

            if x < min(1.0, r)
                # Accept proposed jump
                para_old = para_new
                post_old = post_new
                like_old = like_new
                propdist.μ = para_new

                TTT_old = out["TTT"]
                RRR_old = out["RRR"]
                CCC_old = out["CCC"]

                zend_old = out["zend"]
                ZZ_old = out["ZZ"]
                DD_old = out["DD"]
                QQ_old = out["QQ"]

                if verbose 
                    println("Iteration $j: accept proposed jump")
                end
            else
                # Reject proposed jump
                block_rejections += 1

                if verbose 
                    println("Iteration $j: reject proposed jump")
                end
            end


            # Save every (n_times)th draw

            if j % n_times == 0
                like_sim[j/n_times] = like_old
                post_sim[j/n_times] = post_old
                para_sim[j/n_times, :] = para_old'
                TTT_sim[j/n_times, :] = vec(TTT_old)'
                RRR_sim[j/n_times, :] = vec(RRR_old)'
                CCC_sim[j/n_times, :] = vec(CCC_old)'
                z_sim[j/n_times, :] = vec(zend_old)'
            end
        end

        all_rejections += block_rejections
        block_rejection_rate = block_rejections/(n_sim*n_times)
        println("Block $i rejection rate: $block_rejection_rate")


        ## Once every iblock times, write parameters to a file

        # Calculate starting and ending indices for this block (corresponds to a new chunk in memory)
        block_start = n_sim*(i-n_burn-1)+1
        block_end   = block_start+n_sim-1


        # Write data to file if we're past n_burn blocks
        if i > n_burn
            parasim[block_start:block_end, :] = convert(Matrix{Float32}, para_sim)
            postsim[block_start:block_end, :] = convert(Vector{Float32}, post_sim)
            # likesim[block_start:block_end, :] = convert(Vector{Float32}, like_sim)
            TTTsim[block_start:block_end,:]  = convert(Matrix{Float32}, TTT_sim)
            RRRsim[block_start:block_end,:]  = convert(Matrix{Float32}, RRR_sim)
            zsim[block_start:block_end,:]  = convert(Matrix{Float32}, z_sim)
        end
    end # of block

    close(simfile)

    rejection_rate = all_rejections/(n_blocks*n_sim*n_times)
    println("Overall rejection rate: $rejection_rate")
end # of loop over blocks
