function correction(m::AbstractModel, s0::Array{Float64}, ε0::Array{Float64}, yy::Array{Float64,Float64},P0)
    # s0 is 8xnumParticles
    # ε0 is 3xnumParticles
    # yy is data matrix
    # P0 is what?

    # Compute system and store model parameters
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

    # In Matlab code, they use P0 instead of ε0. What is this?
    s_up = repmat(s0,1,numParticles) + chol(P0)'*randn(numStates,numParticles)

    for t=1:T
        yt = yy[:,t]

        ############ First tempering step / Initialization ################
        ε = randn(numErrors, numParticles)
        s_fore = Φ*s_up + sqrtS2*ε

        # Error for each particle
        perror = repmat(yt-A,1,numParticles)-B*s_fore
        
        # Solve for initial tempering parameter ϕ_1
        init_Ineff_func(φ) = Ineff(φ, 0, yt, perror, H, initialize=1)
        φ_1 = fzero(init_Ineff_func, 0.0000001, 1)
        
        # Calculate initial weights
        density_arr=zeros(numParticles)
        for n=1:numParticles
            density_arr[n]=density(φ_1, 0, yt, perror[:,n], H, initialize=1)
        end
        # Normalize weights
        weights = (density_arr.*weights)./sum(density_arr.*weights)
        # Resampling
        id = multinomial_resampling(weights)
        s_up = s_up[:,id]
        ε_up = ε[:,id]  
        # Reset weights to ones
        weights = ones(numParticles)
        # Initialize likelihood
        lik[t] = lik[t]+log(mean(density_arr.*weights))
       
        # Tempering Initialization
        count = 2
        ϕ_old = ϕ_1

        s_fore = Φ*s_up + sqrtS2*ε_up
        perror = repmat(yt-A, 1, numParticles) - B*s_fore        
        check_InEff=InEff(1.0, φ_1, yt, perror, H)         

        ######### Main Tempering Loop #########
        while (check_InEff>rstar)
            # Increment count
            count = count + 1
            # Define Inefficiency function
            InEff_func(φ) = Ineff(φ, φ_old, yt, perror,H)-rstar
            φ_interval = [φ_old 1.0]
            fphi_interval = [InEff_func(φ_old) InEff_func(1.0)]
            
            # Check to ensure solution exists within interval
            if prod(sign(fphi_interval))==-1
 
                # Set ϕ_new to the solution of the inefficiency function over interval
                φ_new = fzero(Ineff_func,φ_interval)
                check_InEff = InEff(1.0,φ_new,yt,perror,H,initialize=0)
                
                # Update densities given new ϕ
                for i=1:numParticles
                    density_arr[i]= density(φ_new, φ_old, yt, perror[:,i], H, initialize=1) 
                end
                # Normalize weights
                weights = (density_arr.*weights)./(sum(density_arr.*weights))
               
                # Resample particles
                id = multinomial_resampling(weights)
                s_up = s_up[:,id]
                ε_up = ε_up[:,id]
               
                # Set weights back to 1
                weights = ones(numParticles)
                # Update likelihood
                lik[t] = lik[t]+log(mean(density_arr.*weights))
                # Update value for c
                c = c*(0.95 + 0.10*exp(16*(acpt-trgt))/(1+exp(16*(acpt-trgt))))
                
                μ = mean(ε_up,2)
                cov_s = (1/numParticles)*(ε_up-repmat(μ,1,numParticles))*(ε_up - repmat(μ,1,numParticles))'
                checkp = diagm(diag(chol(cov_s)))

                # Mutation Step
                acptVec=Array{Float64,1}
                for i = 1:numParticles
                    new_s, new_ε, acpt = mutation(m,yt, s_particle, ε_particle, A, B, R, Φ,H, sqrtS2, cov_mat)
                    s_fore[i] = new_s
                    ε_fore[i] = new_ε
                    acptVec[i] = acpt
                end
                # Calculate average accept rate
                avgAcpt = mean(acptVec)
                # Get error for all numParticles particles
                perror = repmat(yy-A, 1,numParticles)-B*s_fore
                ϕ_old = ϕ_new
                phi_length[t]+=1
            # If no solution exists within interval, set to rstar
            else 
                check_InEff = rstar
            end
        end
        # End of tempering steps. Begin final round with φ=1
        ϕ_new = 1.0
        # Calculate densities
        for i=1:numParticles
            density_arr[i]= density(φ_new, φ_old, yt, perror[:,i],H,initalize=0)
        end
        # Normalize weights
        weights = (density_arr.*weights)./(sum(density_arr.*weights))
        # Resample
        id = multinomial_resampling(weights)
        s_up = s_up[:,id]
        ε_up = ε_up[:,id]
        # Reset weights to 1
        weights= ones(numParticles)
        # Update likelihood
        lik[t] = lik[t] + log(mean(density_arr.*weights))
        # Update c
        c = c*(0.95+.1*exp(16*(acpt-trgt))/(1+exp(16*(acpt-trgt))))
        # Calculate covariance of s
        μ=mean(eps_up,2)
        cov_s = (1/numParticles)*(eps_up - repmat(μ,1,numParticles))*(eps_up-repmat(μ,1,numParticles))'
        checkp = diagm(diag(chol(cov_s)))
        acptVec=zeros(numParticles)
        # Final round of mutation
        for i=1:numParticles
             new_s, new_ε, acpt = mutation(m, yt, s_particle, ε_particle, A, B, R, Φ, H, sqrtS2, cov_mat)
              s_fore[i] = new_s
              ε_fore[i] = new_ε
              acptVec[i] = acpt 
        end
    # Store for next time iteration
    acptAvg = mean(acptVec)
    Neff[t] = (numParticles^2)/sum(weights.^2)
    s_up = s_fore
    end
    return Neff, lik, acptAvg
end