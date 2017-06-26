function mutation(m::AbstractModel, yt::Array{Float64,1},s_init::Array{Float64,1}, ε_init::Array{Float64,1}, A::Array, B::Array, R::Array, Φ::Array, H::Array, sqrtS2::Array, cov_mat::Array)
    #= 
    This function runs random walk Metropolis Hastings for ONE particle. The caller should loop through all particles and call the method on each one.
    m: the model 
    yt: 1 x numMeasurements vector at time t for all observed y (GDP, inflation, interest rate, etc.)
    s_init: the starting state
    eps_init: the starting epsilon (state error)
    =#

    # Number of Metropolis-Hastings steps. This will eventually be a setting, but hard-coded for now. Should be low.
    N_MH=1
    # Set path
    path = dirname(@__FILE__)
    ind_s=0
    ind_ε=0
    #Get the meta-paramater c from model settings (I think we want a separate c for this part rather than what is used for regular MH but can change later)
    c=get_setting(m,:c)
    # Initialize acceptance counter to zero
    acpt=0
    
    for i=1:N_MH
        # Isolate random matrix for testing purposes
        randMat = randn(size(cov_mat,1),1)
        h5open("$path/../../test/reference/mutationRandomMatrix.h5","w") do file
            write(file, "randMat", randMat)
        end
        # Generate new draw of ε from a N(ε_init, c²cov_mat) distribution
        ε_new=ε_init + c*cov_mat'*randMat
        # Use the state equation to calculate the corresponding state from that ε 
        s_new_fore = Φ*s_init+sqrtS2*ε_new
        # Use the state equation to calculate the state corresponding to ε_init
        s_init_fore = Φ*s_init+sqrtS2*ε_init

        # Calculate difference between data and expected y from measurement equation and calculated states from above for both the new draw of ε_new and the old ε_init (we do this to calculate probabilities below. Since the error is still drawn from a Normal and everything is still linear, we know that y will also be normal. See equation of multivariate normal for how error_new and error_init enter into the pdf).
        error_new = yt-B*s_new_fore-A
        error_init = yt-B*s_init_fore-A
        # Calculate the top and bottom probabilities for the α ratio.
        post_new = log(pdf(MvNormal(zeros(length(yt)),H),error_new)[1]*pdf(MvNormal(zeros(length(ε_new)),eye(length(ε_new),length(ε_new))),ε_new)[1])
        post_init = log(pdf(MvNormal(zeros(length(yt)),H),error_init)[1]*pdf(MvNormal(zeros(length(ε_init)),eye(length(ε_init),length(ε_init))),ε_init)[1])
       
        # α represents the probability of accepting the new particle (post_new and post_init are in logs so subtract and take the exponential
        α = exp(post_new - post_init)
        
        # Accept the particle with probability α (rand() generates Unif(0,1) r.v. If accept set s_init to the particle and ε_init to the error for starting the loop over again
        if rand()<α 
            # Accept
            ind_s = s_new_fore
            ind_ε = ε_new
            acpt = acpt+1
        else 
            # Reject and keep the old particle unchanged
            ind_s = s_init_fore
            ind_ε = ε_init
            #acpt = 0
        end
    end
    return ind_s, ind_ε, acpt 
end

