# This program produces and saves draws from the posterior distribution of
# the parameters.
function estimate{T<:AbstractDSGEModel}(Model::Type{T})
    ### Step 1: Initialization
    model = Model()

    # TODO: load data
    # YY = loaddata()

    post = posterior(model, YY)



    ### Step 2: Finding the posterior mode: re-maximization

    #=
    # Inputs to minimization algorithm (csminwel)
    x0 = invtrans(params, trspec)
    H0 = eye(npara)*1E-4
    nit = 1000
    crit= 1E-10
    MIN = 1
    =#



    ### Step 3: Compute the inverse Hessian
    # Once we find and save the posterior mode, we calculate the Hessian
    # at the mode. The hessian is used to calculate the variance of the proposal
    # distribution, which is used to draw a new parameter in each iteration of
    # the algorithm.

    # Step 3a: Calculate Hessian corresponding to the posterior mode

    # Step 3b: Calculate inverse of Hessian (by Singular Value Decomposition)



    ### Step 4: Initialize algorithm by drawing para_old from a normal distribution centered on the posterior mode (propdens).



    ### Step 5: For nsim*ntimes iterations within each block, generate a new parameter draw. Decide to accept or reject, and save every ntimes_th draw that is accepted.



    ### Step 6: Calculate parameter covariance matrix

end
