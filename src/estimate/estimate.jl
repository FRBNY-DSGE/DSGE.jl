# This program produces and saves draws from the posterior distribution of
# the parameters.
function estimate{T<:AbstractModel}(Model::Type{T})
    ### Step 1: Initialize model
    model = Model()
    spec = model.spec

    # TODO: load data
    mf = MatFile("$savepath/YY.mat")
    YY = get_variable(mf, "YY")
    close(mf)

    post = posterior(model, YY)

    
    
    ### Step 2: Find posterior mode

    # Specify starting mode
    println("Reading in previous mode")
    mf = MatFile("$savepath/mode_in.mat")
    mode = get_variable(mf, "params")
    close(mf)
    update!(model.Θ, mode)

    if spec["reoptimize"]
        println("Reoptimizing")
        
        # Inputs to minimization algorithm
        function posterior_min!{T<:FloatingPoint}(x::Vector{T})
            tomodel!(x, model.Θ)
            return -posterior(model, YY)
        end

        xh = toreal(model.Θ)
        H = 1e-4 * eye(spec["n_params"])
        nit = 1000
        crit = 1e-10    
        converged = false

        # If the algorithm stops only because we have exceeded the maximum number of
        # iterations, continue improving guess of modal parameters
        while !converged
            out, H = csminwel(posterior_min!, xh, H; ftol=crit, iterations=nit, show_trace=true, verbose=true)
            xh = out.minimum
            converged = !out.iteration_converged
        end

        # Transform modal parameters so they are no longer bounded (i.e., allowed
        # to lie anywhere on the real line).
        tomodel!(xh, model.Θ)
        mode = [α.value for α in model.Θ]

        # Write mode to file
        mf = MatFile("$savepath/mode_out.mat", "w")
        put_variable(mf, "mode", mode)
        close(mf)
    end    
    

    
    ### Step 3: Compute proposal distribution

    # Calculate the Hessian at the posterior mode
    if spec["recalculate_hessian"]
        println("Recalculating Hessian")
        hessian, _ = hessizero!(mode, model, YY; noisy=true)

        mf = MatFile("$savepath/hessian.mat", "w")
        put_variable(mf, "hessian", hessian)
        close(mf)
    else
        println("Using pre-calculated Hessian")
        mf = MatFile("$savepath/hessian.mat")
        hessian = get_variable(mf, "hessian")
        close(mf)
    end
    
    # The hessian is used to calculate the variance of the proposal
    # distribution, which is used to draw a new parameter in each iteration of
    # the algorithm.
    propdist = proposal_distribution(mode, hessian)
    if propdist.rank != spec["n_free_params"]
        println("problem – shutting down dimensions")
    end

    

    ### Step 4: Metropolis-Hastings algorithm

    # Set the jump size for sampling
    cc0 = 0.01
    cc = 0.09
    
    metropolis_hastings(propdist, model, YY, cc0, cc)


    
    ### Step 5: Calculate parameter covariance matrix
    
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

    return DegenerateMvNormal(μ, σ)
end



function metropolis_hastings{T<:FloatingPoint}(propdist::Distribution, model::AbstractModel, YY::Matrix{T}, cc0::T, cc::T, randvecs = [], randvals = [])

    # If testing, then we read in a specific sequence of "random" vectors and numbers
    testing = !(randvecs == [] && randvals == [])
    println("Testing = $testing")
    
    # Initialize algorithm by drawing para_old from a normal distribution centered on the posterior mode until the parameters are within bounds or the posterior value is sufficiently large    
    para_old = propdist.μ
    post_old = -Inf
    initialized = false
    while !initialized
        if testing
            para_old = propdist.μ + cc0*propdist.σ*randvecs[:, 1]
        else
            para_old = rand(propdist; cc=cc0)
        end
        
        post_old, like_old, out = posterior!(para_old, model, YY; mh=true)
        if post_old > -Inf
            propdist.μ = para_old
            initialized = true
        end
    end

    # For n_sim*n_times iterations within each block, generate a new parameter draw. Decide to accept or reject, and save every (n_times)th draw that is accepted.
    spec = model.spec

    if testing
        n_blocks = 1
        n_sim = 200
        n_times = 5
    else
        n_blocks = spec["n_blocks"]
        n_sim = spec["n_sim"]
        n_times = spec["n_times"]
    end

    all_rejections = 0

    para_sim = zeros(n_sim, spec["n_params"])
    like_sim = zeros(n_sim)
    post_sim = zeros(n_sim)
    TTT_sim = zeros(n_sim, spec["n_states_aug"]^2)
    RRR_sim = zeros(n_sim, spec["n_states_aug"]*spec["n_exoshocks"])
    CCC_sim = zeros(n_sim, spec["n_states_aug"])
    z_sim = zeros(n_sim, spec["n_states_aug"])
    
    for i = 1:n_blocks
        block_rejections = 0

        for j = 1:(n_sim*n_times)
            # Draw para_new from the proposal distribution
            if testing
                k = (i-1)*(n_sim*n_times) + j
                para_new = propdist.μ + cc*propdist.σ*randvecs[:, k]
            else
                para_new = rand(propdist; cc=cc)
            end
            
            # Solve the model, check that parameters are within bounds, and
            # evaluate the posterior.
            post_new, like_new, out = posterior!(para_new, model, YY; mh=true)
            if testing
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
                x = randvals[k]
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

                if testing
                    println("Iteration $j: accept proposed jump")
                end
            else
                # Reject proposed jump
                block_rejections += 1
                
                if testing
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
    end
    
    rejection_rate = all_rejections/(n_blocks*n_sim*n_times)
    println("Overall rejection rate: $rejection_rate")
end
