# This program produces and saves draws from the posterior distribution of
# the parameters.
function estimate{T<:AbstractModel}(Model::Type{T})
    ### Step 1: Initialize model
    model = Model990()
    spec = model.spec

    # TODO: load data
    # YY = loaddata()

    post = posterior(model, YY)

    
    
    ### Step 2: Find posterior mode

    # Specify starting mode
    mode_in = readcsv("$savepath/mode_in.csv")
    update!(model, mode_in)

    # Inputs to minimization algorithm (csminwel)
    x0 = [toreal(α) for α in model.Θ]
    H0 = 1e-4 * eye(spec["n_params"])
    nit = 1000
    crit = 1e-10

    # We seek the posterior distribution of the parameters.
    # We first find the mode of the distribution (i.e., the maximum of its pdf)
    # so that once we begin sampling, we can begin from this mode.
    xh = csminwel(posterior!, x0, H0, model, YY; ftol=crit, iterations=nit)

    # Transform modal parameters so they are no longer bounded (i.e., allowed
    # to lie anywhere on the real line).
    update!(model, xh)
    mode = [tomodel(α) for α in model.Θ]
    update!(model, mode)
    writecsv("$savepath/mode_out.csv", mode)

    # If the algorithm stops only because we have exceeded the maximum number of
    # iterations, continue improving guess of modal parameters by recursively
    # calling gibb.
    
    

    ### Step 3: Compute proposal distribution

    # Calculate the Hessian at the posterior mode
    hessian, _ = hessizero!(mode, model, YY; noisy=true)
    writecsv("$savepath/hessian.csv", hessian)
    
    # Set the jump size for sampling
    cc0 = 0.01
    cc = 0.09
    date_q = 1:spec["quarters_ahead"]

    # The hessian is used to calculate the variance of the proposal
    # distribution, which is used to draw a new parameter in each iteration of
    # the algorithm.
    propdist = proposal_distribution(mode, hessian, cc0)
    if propdist.rank != spec["n_free_params"]
        println("problem – shutting down dimensions")
    end

    

    ### Step 4: Metropolis-Hastings algorithm



    ### Step 5: Calculate parameter covariance matrix
            
end



# Compute proposal distribution: degenerate normal with mean μ and covariance cc*hessian^(-1)
function proposal_distribution{T<:FloatingPoint}(μ::Vector{T}, hessian::Matrix{T}, cc::T)
    n = size(hessian, 1)
    @assert (n, n) = size(hessian)

    u, s_diag, v = svd(hessian)
    Σ_inv = hessian/cc^2
    big_evals = find(x -> x > 1e-6, s_diag)
    rank = length(big_evals)

    logdet = 0.0
    Σ = zeros(n, n)
    for i = 1:rank
        Σ[i, i] = cc^2/s_diag[i]
        logdet += log(Σ[i, i])
    end

    σ = cc*u*sqrt(Σ)

    return DegenerateMvNormal(μ, σ, Σ, Σ_inv, rank, logdet)
end



function metropolis_hastings{T<:FloatingPoint}(mode::Vector{T}, propdist::Distribution, model::AbstractModel, YY::Matrix{T})
    # Initialize algorithm by drawing para_old from a normal distribution centered on the posterior mode until the parameters are within bounds or the posterior value is sufficiently large    
    para_old = mode
    initialized = false
    while !initialized
        para_old = rand(propdist)
        post_old = posterior!(para_old, model, YY; checkbounds=true)
        #propdens = logpdf(propdist, para_old)
        if post_old > -Inf
            initialized = true
        end
    end

    # Initialize variables for calculating rejection rate
    Tim = 0
    eT = 0
    reje = 0

    # For n_draws*ntimes iterations within each block, generate a new parameter draw. Decide to accept or reject, and save every ntimes_th draw that is accepted.
    for i = 1:spec["n_blocks"]
        parasim = zeros(spec["n_draws"], spec["n_params"])
        likesim = zeros(spec["n_draws"], 1)
        postsim = zeros(spec["n_draws"], 1)
        rej = zeros(spec["n_draws"]*spec["n_times"], 1)
        TTTsim = zeros(spec["n_draws"], spec["n_states_aug"]^2)
        RRRsim = zeros(spec["n_draws"], spec["n_states_aug"]*spec["n_exoshocks"])
        CCCsim = zeros(spec["n_draws"], spec["n_states_aug"])
        zsim = zeros(spec["n_draws"], spec["n_states_aug"])

        for j = 1:(spec["n_draws"]*spec["n_times"])
            Tim += 1

            # Draw para_new from the proposal distribution
            para_new = rand(propdist)
            
            # Solve the model, check that parameters are within bounds, and
            # evaluate the posterior.
            post_new = posterior!(para_new, model, YY; checkbounds=true)

            [post_new,like_new,zend_new,ZZ_new,DD_new,QQ_new] = feval('objfcnmhdsge',para_new,bounds,YY,YY0,nobs, nlags,nvar,mspec,npara,trspec,pmean,pstdd,pshape,TTT_new,RRR_new,CCC_new,valid_new,para_mask,coint,cointadd,cointall,YYcoint0,args_nant_antlags{:});

            # Calculate the multivariate log likelihood of jump from para_old to para_new
            propdens = logpdf(propdist, para_new)
        end
    end
end
