function mutation(m::AbstractModel, yt::Array{Float64,1},s_init::Array{Float64,1}, ε_init::Array{Float64,1}, A::Array, B::Array, R::Array, Φ::Array, H::Array, sqrtS2::Array, cov_mat::Array,N_MH::Int64)
#=This function runs random walk Metropolis Hastings for ONE particle. The caller should loop through all particles and call the method on each one.
m: The model. Used to get the setting tpf_c (which is VERY important it turns out...) 
yt: A 1 x num_measurement vector at time t for all observed y (GDP, PCE, FFR, etc.)
s_init: The starting state before mutation.
ε_init: The starting epsilon (state error) before mutation.
=#

    # Set path--used only for testing
    path = dirname(@__FILE__)
    #Initialize ind_s and ind_ε
    ind_s=s_init
    ind_ε=ε_init
    #Get the meta-paramater tpf_c from model setting. This updates dynamically to tune TPF over time (which is very important as we all know).
    c=get_setting(m,:tpf_c)
    # Initialize acceptance counter to zero
    acpt=0
    
    for i=1:N_MH
        # Isolate random matrix for testing purposes
        rand_mat = randn(size(R,1),1)
        h5open("$path/../../test/reference/mutationRandomMatrix.h5","w") do file
            write(file, "randMat", randMat)
        end
        # Generate new draw of ε from a N(ε_init, c²R) distribution where R is cov_s (defined in tpf.jl) and c is the tuning parameter
        ε_new=ε_init + c*Matrix(chol(nearestSPD(R)))'*rand_mat
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

        # α represents the probability of accepting the new particle (post_new and post_init are in logs so subtract and take the exponential)
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
#THIS IS DIFFERENT THAN MATLAB BUT WE THINK MATLAB WRONG
#            acpt = 0
        end
    end
    acpt /= N_MH
    return ind_s, ind_ε, acpt 
end

