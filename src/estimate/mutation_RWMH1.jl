function mutation_RWMH1(m::AbstractModel, yt::Array{Float64,1},s_init, ε_init)
    #=This function runs random walk metropolis hastings for ONE particle. The caller should loop through all particles and call the method on each one.
    m: the model 
    yt: 1xnumMeasurements vector at time t for all observed y (GDP, inflation, interest rate, etc.)
    s_init: the starting state
    eps_init: the starting epsilon (state error)
    =#

    #number of metropolis-hastings steps. this will eventually be a setting, but hard-coded for now. Should be low.
    N_MH=1
   
    #solve for the matrices in the state equation (using TPF paper notation, not standard)
    φ,R,C=solve(m)
    #solve for matrices in measurement equation
    measure = measurement(m,φ,R,C)
    #H is the covariance matrix of the measurement equation
    H=measure[:EE]
    
    #varStateEq is the covariance matrix of the state equation(should change name. TPF paper calls it S2 but probably want to call it something different.
    varStateEq=measure[:MM]*measure[:QQ]*measure[:MM]'
    #Get the meta-paramater c from model settings (I think we want a separate c for this part rather than what is used for regular MH but can change later)
    c=get_setting(m,:c)
    #Convert the varStateEq matrix into a covariance matrix of the state error
    U,E,V=svd(varStateEq)
    cov_mat=U*diagm(sqrt(E))
    #initialize acceptance counter to zero
    acpt=0
    
    for i=1:N_MH
        #Generate new draw of ε from a N(ε_init, c²cov_mat) distribution
        ε_new=ε_init + c*cov_mat*randn(3)
        #Use the state equation to calculate the corresponding state from that ε 
        s_new_fore = φ*s_init+R*varStateEq*ε_new
        #Use the state equation to calculate the state corresponding to ε_init
        s_init_fore = φ*s_init+R*varStateEq*ε_init

        #Calculate difference between data and expected y from measurement equation and calculated states from above for both the new draw of ε_new and the old ε_init (we do this to calculate probabilities below. Since the error is still drawn from a Normal and everything is still linear, we know that y will also be normal. See equation of multivariate normal for how error_new and error_init enter into the pdf).
        error_new = yt-measure[:ZZ]*s_new_fore-measure[:DD]
        error_init = yt-measure[:ZZ]*s_init_fore - measure[:DD]
        #Calculate the top and bottom probabilities for the α ratio.
        post_new = log(pdf(MvNormal(zeros(length(yt)),H),error_new)*pdf(MvNormal(zeros(length(ε_new)),eye(length(ε_new),length(ε_new))),ε_new))
        post_init = log(pdf(MvNormal(zeros(length(yt)),H),error_init)*pdf(MvNormal(zeros(length(ε_init)),eye(length(ε_init),length(ε_init))),ε_init))
       
        #α represents the probability of accepting the new particle (post_new and post_init are in logs so subtract and take the exponential
        α = exp(post_new - post_init)
        
        #Accept the particle with probability α (rand() generates Unif(0,1) r.v. If accept set s_init to the particle and ε_init to the error for starting the loop over again
        if rand()<α 
            #accept
            s_init = s_new_fore
            ε_init = ε_new
            acpt = acpt+1
        #else reject and keep the old particle unchanged
        end
    end
    return s_init, ε_init, acpt 
end

