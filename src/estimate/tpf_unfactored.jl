using DSGE
using Roots
function tpf(m::AbstractModel, yy::Array, s0::Array{Float64}, P0::Array)
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
    cov_mat = getChol(varStateEq)
    sqrtS2=R*eye(size(cov_mat',1)) # TEMPORARY CHANGE - for TESTING (old code: sqrtS2=R*cov_mat'
    
    # Get tuning parameters from the model
    rstar = get_setting(m,:tpf_rstar)
    c = get_setting(m,:tpf_c)
    acpt_rate = get_setting(m,:tpf_acpt_rate)
    trgt = get_setting(m,:tpf_trgt)
    N_MH = get_setting(m,:tpf_N_MH)
    numParticles = get_setting(m,:tpf_numParticles)

    # Initialize time
    T = size(yy,2)
    # Size of covariance matrix
    numErrors = size(sqrtS2,2)
    #number of states
    numStates = size(B, 2)
    
    # Likelihood matrix
    lik = zeros(T)
    Neff=zeros(T)
    len_phis = ones(T)
    weights = ones(numParticles)
    density_arr=zeros(numParticles)

    # In Matlab code, they use P0 instead of ε0. What is this?
    s_up = repmat(s0,1,numParticles) + getChol(P0)'*randn(numStates,numParticles)

    for t=1:T
        yt = yy[:,t]

        ############ First tempering step / Initialization ################
        ε = randn(numErrors, numParticles)
        s_fore = Φ*s_up + sqrtS2*ε

        # Error for each particle
        perror = repmat(yt-A,1,numParticles)-B*s_fore
        # Solve for initial tempering parameter ϕ_1
        init_Ineff_func(φ) = ineff_func(φ, 0.0, yt, perror, H, initialize=1)-rstar
        φ_1 = fzero(init_Ineff_func,0.000001, 1.0)
        
        # Update weights array and resample particles.
        loglik, weights, s_up, ε_up = correctAndResample(φ_1, 0.0,yt,perror,density_arr,weights,s_up,ε,H,numParticles,initialize=1)
        # Update likelihood
        lik[t]=lik[t]+loglik

        #= # Calculate initial weights
        density_arr=zeros(numParticles)
        for n=1:numParticles
            density_arr[n]=density(φ_1, 0.0, yt, perror[:,n], H, initialize=1)
            if isnan(density_arr[n])
                println("NaN in density array!")
            end
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
=#
        # Tempering Initialization
        count = 2
        φ_old = φ_1

        s_fore = Φ*s_up + sqrtS2*ε_up
        perror = repmat(yt-A, 1, numParticles) - B*s_fore        
        check_InEff=ineff_func(1.0, φ_1, yt, perror, H)         

        ######### Main Tempering Loop #########
        while (check_InEff>rstar)
            # Increment count
            count = count + 1
            # Define Inefficiency function
            InEff_func(φ) = ineff_func(φ, φ_old, yt, perror,H)-rstar
            φ_interval = [φ_old, 1.0]
            fphi_interval = [InEff_func(φ_old) InEff_func(1.0)]
            
            # Check to ensure solution exists within interval
            if prod(sign(fphi_interval))==-1
               
                # Set ϕ_new to the solution of the inefficiency function over interval
                φ_new = fzero(InEff_func,φ_interval)
                check_InEff = ineff_func(1.0,φ_new,yt,perror,H,initialize=0)
                # Update weights array and resample particles.
                loglik, weights, s_up, ε_up = correctAndResample(φ_new,φ_old,yt,perror,density_arr,weights,s_up,ε_up,H,numParticles,initialize=0)
                # Update likelihood
                lik[t]=lik[t]+loglik

               #= 
                # Update densities given new ϕ
                for i=1:numParticles
                    density_arr[i]= density(φ_new, φ_old, yt, perror[:,i], H, initialize=0)
                end
                # Normalize weights
                weights = (density_arr.*weights)./(sum(density_arr.*weights))
                println("Within while-loop weights: ")
                println(weights)
                # Resample particles
                id = multinomial_resampling(weights)
                s_up = s_up[:,id]
                ε_up = ε_up[:,id]
               
                # Set weights back to 1
                weights = ones(numParticles)
                # Update likelihood
                lik[t] = lik[t]+log(mean(density_arr.*weights))
                =#

                # Update value for c
                c = updateC(c,acpt_rate,trgt)
                
                μ = mean(ε_up,2)
                cov_s = (1/numParticles)*(ε_up-repmat(μ,1,numParticles))*(ε_up - repmat(μ,1,numParticles))'
                checkp = diagm(diag(getChol(cov_s)))

                # Mutation Step
                acptVec=zeros(numParticles)
                
                for i = 1:numParticles  
                    new_s, new_ε, acpt = mutation(m,yt,s_up[:,i],ε_up[:,i],A,B,cov_s,Φ,H,sqrtS2,cov_mat,N_MH)
                    s_fore[:,i] = new_s
                    ε_up[:,i] = new_ε
                    acptVec[i] = acpt
                end
                # Calculate average accept rate
                acpt_rate = mean(acptVec)
                # Get error for all numParticles particles
                perror = repmat(yt-A, 1,numParticles)-B*s_fore
                φ_old = φ_new
                len_phis[t]+=1
            # If no solution exists within interval, set to rstar
            else 
                check_InEff = rstar
            end
        end
        # End of tempering steps. Begin final round with φ=1
        φ_new = 1.0
        # Calculate densities
        for i=1:numParticles
            density_arr[i]= density(φ_new, φ_old, yt, perror[:,i],H,initialize=0)
            if density_arr[i]==0.0
                println("Density array is zero at this index!")
                println("φ_new")
                println(φ_new)
                println("φ_old")
                println(φ_old)
                println("perror[:,i]")
                println(perror[:,i])
                println("inv(H)")
                println(inv(H))
             end
        end
        # ERROR ARISES HERE - perror gets too large, causing density to go to 0.
        println("Last step density arr")
        println(density_arr)
        println("last step old weights")
        println(weights)

        # Normalize weights
        weights = (density_arr.*weights)./(sum(density_arr.*weights))
        
        # After ~60th iteration, weights -> NaN
        println("Last step weights vector")
        println(weights)

        # Resample
        id = multinomial_resampling(weights)
        s_up = s_up[:,id]
        ε_up = ε_up[:,id]
        # Reset weights to 1
        weights = ones(numParticles)
        # Update likelihood
        lik[t] = lik[t] + log(mean(density_arr.*weights))
        # Update c
        c = updateC(c,acpt_rate,trgt)
        # Calculate covariance of s
        μ=mean(ε_up,2)
        cov_s = (1/numParticles)*(ε_up - repmat(μ,1,numParticles))*(ε_up-repmat(μ,1,numParticles))'
        checkp = diagm(diag(getChol(cov_s)))
        acptVec=zeros(numParticles)
        # Final round of mutation
        for i=1:numParticles
            new_s, new_ε, acpt = mutation(m, yt, s_up[:,i], ε_up[:,i], A, B, R, Φ, H, sqrtS2, cov_mat,N_MH)
            s_fore[:,i] = new_s
            ε_up[:,i] = new_ε
            acptVec[i] = acpt 
        end
        # Store for next time iteration
        acpt_rate = mean(acptVec)
        Neff[t] = (numParticles^2)/sum(weights.^2)
        s_up = s_fore
    println(t)
    end
    return Neff, lik
end



# Function for calculating cholesky decomposition of matrix
function getChol(mat::Array)
    U,E,V=svd(mat)
    return U*diagm(sqrt(E))
end

# Calculation for updating c
function updateC(c::Float64, acpt::Float64, trgt::Float64)
    return c*(0.95+0.1*exp(16*(acpt-trgt))/(1+exp(16*(acpt-trgt))))
end

################
function correctAndResample(φ_new, φ_old, yt, perror,density_arr,weights, s_up, ε, H, numParticles; initialize=0)
    # Calculate initial weights
    for n=1:numParticles
        density_arr[n]=density(φ_new, φ_old, yt, perror[:,n], H, initialize=initialize)
    end
    # Normalize weights
    weights = (density_arr.*weights)./sum(density_arr.*weights)
    # Resampling
    id = multinomial_resampling(weights)
    s_up = s_up[:,id]
    ε_up = ε[:,id]  
    # Reset weights to ones
    weights = ones(numParticles)
    # Calculate likelihood
    loglik = log(mean(density_arr.*weights))

    return loglik, weights, s_up, ε_up
end
##############3
