# This program produces and saves draws from the posterior distribution of
# the parameters.
function estimate{T<:AbstractModel}(Model::Type{T})
    ### Step 1: Initialization
    model = Model990()
    spec = model.spec

    # TODO: load data
    # YY = loaddata()

    post = posterior(model, YY)

    
    
    ### Step 2: Finding the posterior mode: re-maximization

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

    # Once we get a mode, we save the parameter vector
    writecsv("$savepath/mode_out.csv", mode)

    # If the algorithm stops only because we have exceeded the maximum number of
    # iterations, continue improving guess of modal parameters by recursively
    # calling gibb.
    
    

    ### Step 3: Compute the inverse Hessian
    # Once we find and save the posterior mode, we calculate the Hessian
    # at the mode. The hessian is used to calculate the variance of the proposal
    # distribution, which is used to draw a new parameter in each iteration of
    # the algorithm.

    ## Step 3a: Calculate Hessian corresponding to the posterior mode
    hessian, _ = hessizero!(mode, model, YY; noisy=true)

    # Save computed Hessian in output file
    writecsv("$savepath/hessian.csv", hessian)

    
    ## Step 3b: Calculate inverse of Hessian (by Singular Value Decomposition)
    rankhhm = spec["n_params"]
    u, s_diag, v = svd(hessian)
    sigpropinv = copy(hessian)
    big_evals = find(x -> x > 1e-6, s_diag)
    sigpropdim = length(big_evals)

    sigproplndet = 0.0
    s = zeros(spec["n_params"], spec["n_params"])
    for i = 1:spec["n_params"]
        if i > sigpropdim
            s[i, i] = 0.0
        else
            s[i,i] = 1/s_diag[i]
            sigproplndet += log(s[i, i])
        end
    end

    invhhm  = u*s*u'
    sigscale = u*sqrt(s)

    if length(bigev) != spec["n_free_params"]
        println("problem – shutting down dimensions")
    end

    

    ### Step 4: Initialize algorithm by drawing para_old from a normal distribution centered on the posterior mode (propdens).



    ### Step 5: For nsim*ntimes iterations within each block, generate a new parameter draw. Decide to accept or reject, and save every ntimes_th draw that is accepted.



    ### Step 6: Calculate parameter covariance matrix
            
end
