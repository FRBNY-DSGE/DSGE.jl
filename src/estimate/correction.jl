function correction(m::AbstractModel, s0::Array{Float64}, ε0::Array{Float64}, yy::Array{Float64,Float64})
    #s is 8xnumParticles
    #ε is 3xnumParticles
    #φn is scalar
    #φt1 is scalar (phi n-1)
    #w is a vector of weights
    #B is psi, the measurement error function (here since linear, we just premultiply it by s
    
    #s_up is what?

    #Get all the stuff of import from the model
    sys=compute_system(m)
    R=sys.transition.RRR
    Φ=sys.transition.TTT
    C=sys.transition.CCC
    H=sys.measurement.EE
    A=sys.measurement.DD
    B=sys.measurement.ZZ
    varStateEq=sys.measurement.MM*sys.measurement.QQ*sys.measurement.MM'
    U,E,V=svd(varStateEq)
    cov_mat = U*diagm(sqrt(E))
    sqrtS2=R*cov_mat'

    # Initialize time
    T = size(yy,2)
    # Size of covariance matrix
    numErrors = size(S2,1)
    #number of states
    numStates = size(B, 2)
    # Number of particles
    numParticles= length(s)
    
    # Likelihood matrix
    lik = zeros(T)
    Neff=zeros(T)
    len_phis = ones(T)
    weights = ones(numParticles)

    #In Matlab code, use P0 instead of ε0. what is this?
    s_up = repmat(s0,1,numParticles) + chol(ε0)'*randn(numStates,numParticles)

   
    for t=1:T
        yt = yy[:,t]

        #first temper step
        ε = randn(numErrors, numParticles)
        s_fore = Φ*s_up + sqrtS2*ε

        #error for each particle
        perror=repmat(yt-A,1,N)-B*s_fore
        
        # Initial tempering step with φ=φ_initial
        init_Ineff_func(φ) = Ineff(φn, φn_1, yt, perror, H, 1)
        φ_1=fzero(init_Ineff_func, 0.0000001, 1)
        
        # Calculate initial weights
        density_arr=zeros(numParticles)
        for n=1:numParticles
            density_arr[n]=density(φ_1, 0, yt, perror[:,n],H,initialize=1)
        end

        weights = (density_arr.*weights)./sum(density_arr.*weights)

        id = multinomial_resampling(weights)
        s_up = s_up[:,id]
        ε_up = ε[:,id]
        
        weights = ones(numParticles)
        # Initialize likelihood
        lik[t] = lik[t]+log(mean(density_arr.*weights))

        # Tempering Initialization
        count = 2
        check_InEff=InEff(1.0, φ1, yt, perror, H) 
        s_fore = Φ*s_up + sqrtS2*ε_up 
        perror = repmat(yt-A, 1, M) - B*s_fore        
        φ_old = φ1
        
        # Tempering Loop
        while (check_InEff>rstar)
            count = count + 1
            InEff_func(φ) = Ineff(φ, φ_old, yt, perror,H)-rstar
            
            φ_interval = [φ_old 1.0]
            
            fphi_interval = [InEff_func(φ_old) InEff_func(1.0)]
            
            if prod(sign(fphi_interval))==-1
                φ_new = fzero(Ineff_func,φ_interval)
                check_InEff = InEff(1.0,φ_new,yt,perror,H, initialize=0)
                
                #s_fore = Φ*s + sqrtS2*ε 
                #perror = repmat(yt-A, 1, M)-B*s
                
                for i=1:numParticles
                    density_arr[i]= density(φn, φn_1, yt, perror[:,i], H, initialize=1) 
                end
            
                weights = (density_arr.*w)./(sum(density_arr.*weights))
                
                id = multinomial_resampling(weights)
                s_up = s_up[:,id]
                ε_up = ε_up[:,id]
                
                # We need to just have this in a for-loop... otherwise space-inefficient 
                weights = ones(numParticles,1)
                lik[t] = lik[t]+log(mean(density_arr.*weights))
            
                c = c*(0.95 + 0.10*exp(16*(acpt-trgt))/(1+exp(16*(acpt-trgt))))
                #This syntax may be incorrect - either for getting the mean or what eps_up translates to
                μ = mean(ε_up,2)
                cov_s = (1/N)*(ε_up-repmat(μ,1,N))*(ε_up - repmat(μ,1,N))'
                checkp = diagm(diag(chol(A)))
                # WE ARE FEELING GOOD ABOUT THE PROSPECTS OF CHECKP BEING UPPER TRIANGULAR
                # if checkp ~= 0
                #    R = diag(diag(cov_s))
                # else
                #    R = cov_s

                #mutation
                acptVec=Array{Float64,1}
                for i = 1:numParticles
                    new_s, new_ε, acpt = mutation(m,yt, s_particle, ε_particle, A, B, R, Φ,H, sqrtS2, cov_mat)
                    s_fore[i] = new_s
                    ε_fore[i] = new_ε
                    acptVec[i] = acpt
                end
                avgAcpt = mean(acptVec)
                #get error for all numParticles particles
                perror = repmat(yy-A, 1,numParticles)-B*s_fore
                phi_old = phi_new
                phi_length[t]+=1
                print("φ=$phi_old")
            else 
                check_InEff = rstar
            end
        end
        #finish of tempering. Begin final round with φ=1
        phi_new = 1.0
        for i=1:numParticles
            density_arr[i]= density(φn, φn_1, yt, perror[:,i],H,initalize=0)
        end
            
        weights = (density_arr.*w)./(sum(density_arr.*weights))
        
        id = multinomial_resampling(weights)
        s_up = s_up[:,id]
        ε_up = ε_up[:,id]
        weights= ones(numParticles)
        lik[t] = lik[t] + log(mean(density_arr.*weights))
        c = c*(0.95+.1*exp(16*(acpt-trgt))/(1+exp(16*(acpt-trgt))))
        μ=mean(eps_up,2)
        cov_s = (1/numParticles)*(eps_up - repmat(μ,1,numParticles))*(eps_up-repmat(μ,1,numParticles))'
        checkp = diagm(diag(chol(A)))
        acptVec=zeros(numParticles)
        for i=1:numParticles
             new_s, new_ε, acpt = mutation(m,yt, s_particle, ε_particle, A, B, R, Φ, H, sqrtS2, cov_mat)
              s_fore[i] = new_s
              ε_fore[i] = new_ε
              acptVec[i] = acpt 
        end
    acptAvg = mean(acptVec)
    Neff[t] = (numParticles^2)/sum(weights.^2)
    s_up = s_fore
    end
    return Neff, lik
end