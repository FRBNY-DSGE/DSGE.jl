using DSGE
using Roots
function tpf(m::AbstractModel, yy::Array, s0::Array{Float64}, P0::Array, A, B, H, R, S2, Φ)
    # s0 is 8xnumParticles
    # ε0 is 3xnumParticles
    # yy is data matrix
    # P0 is what?

####THIS NEEDS TO COME BACK. TAKEN OUT SO THAT WE CAN PASS IN MATRICES FROM SHORFHEIDGE MATLAB CODE, BUT WILL NEED TO CHANGE
    #Compute system and store model parameters
  #=  sys=compute_system(m)
    R=sys.transition.RRR
    Φ=sys.transition.TTT
    C=sys.transition.CCC
    H=sys.measurement.EE
    A=sys.measurement.DD
    B=sys.measurement.ZZ
    S2=sys.measurement.QQ
    =#   
    sqrtS2=R*getChol(S2)'

    # Get tuning parameters from the model
  
    ###WILL ALSO NEED TO BRING THIS BACK LATER
    #= rstar = get_setting(m,:tpf_rstar)
    c = get_setting(m,:tpf_c)
    acpt_rate = get_setting(m,:tpf_acpt_rate)
    trgt = get_setting(m,:tpf_trgt)
    N_MH = get_setting(m,:tpf_N_MH)
    numParticles = get_setting(m,:tpf_numParticles)
    =#
    rstar= 2.0
    acpt_rate= 0.5
    trgt= 0.25
    c = 0.1
    N_MH= 2
    numParticles=100


    # End time (last period)
    T = size(yy,2)
    # Size of covariance matrix
    numErrors = size(S2,1)
    #number of states
    numStates = size(B, 2)
    
    # Likelihood matrix
    lik = zeros(T)
    Neff=zeros(T)
    len_phis = ones(T)
    weights = ones(numParticles)
    density_arr=zeros(numParticles)
    ε_up = zeros(numErrors)

    # In Matlab code, they use P0 instead of ε0. What is this?
    s_up = repmat(s0,1,numParticles) + getChol(P0)'*randn(numStates,numParticles)
    
    acptCheck = ones(T)
    alph=ones(numParticles)

    for t=1:T
        yt = yy[:,t]

        ############ First tempering step / Initialization ################
        ε = randn(numErrors, numParticles)
        s_fore = Φ*s_up + sqrtS2*ε

        # Error for each particle
        perror = repmat(yt-A,1,numParticles)-B*s_fore

        ### Solve for initial tempering parameter ϕ_1
        
        #function to solve for zeros
        init_Ineff_func(φ) = ineff_func(φ, 2.0*pi, yt, perror, H, initialize=1)-rstar
        #solve for zeros
        φ_1 = fzero(init_Ineff_func,0.000001, 1.0)

        # Update weights array and resample particles.
        loglik, weights, s_up, ε_up = correctAndResample(φ_1,0.0,yt,perror,density_arr,weights,s_up,ε,H,numParticles,initialize=1)
        # Update likelihood
        lik[t]=lik[t]+loglik

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
                @show acpt_rate
                # Update value for c
                c = updateC(m,c,acpt_rate,trgt)
 
                @show c
                μ = mean(ε_up,2)
                cov_s = (1/numParticles)*(ε_up-repmat(μ,1,numParticles))*(ε_up - repmat(μ,1,numParticles))'
 
                # Mutation Step
                acptVec=zeros(numParticles)
                
                for i = 1:numParticles
                    ind_s, ind_ε, ind_acpt,α = mutation(m,yt,s_up[:,i],ε_up[:,i],A,B,cov_s,Φ,H,sqrtS2,S2,N_MH)
                    s_fore[:,i] = ind_s
                    ε_up[:,i] = ind_ε
                    acptVec[i] = ind_acpt
                    alph[i]=α
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
        # Update weights array and resample particles.
        loglik, weights, s_up, ε_up = correctAndResample(φ_new,φ_old,yt,perror,density_arr,weights,s_up,ε_up,H,numParticles,initialize=0)
        # Update likelihood
        lik[t]=lik[t]+loglik
        
        c = updateC(m,c,acpt_rate,trgt)

        μ = mean(ε_up,2)
        cov_s = (1/numParticles)*(ε_up-repmat(μ,1,numParticles))*(ε_up - repmat(μ,1,numParticles))'
#=
        # Calculate densities
        for i=1:numParticles
            density_arr[i]= density(φ_new, φ_old, yt, perror[:,i],H,initialize=0)
            if density_arr[i]==0.0
                printlmn("Density array is zero at this index!")
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
=#
        acptVec=zeros(numParticles)
        # Final round of mutation
        for i=1:numParticles
            ind_s, ind_ε, ind_acpt,α = mutation(m, yt, s_up[:,i], ε_up[:,i], A, B, cov_s, Φ, H, sqrtS2,S2,N_MH)
            s_fore[:,i] = ind_s
            ε_up[:,i] = ind_ε
            acptVec[i] = ind_acpt 
            alph[i] = α
        end
        # Store for next time iteration
        acpt_rate = mean(acptVec)
        acptCheck[t]=acpt_rate
        Neff[t] = (numParticles^2)/sum(weights.^2)
        s_up = s_fore
    println("Period: $t")
    end

    @show mean(acptCheck)
    @show  mean(alph)
    return Neff, lik
end



# Function for calculating cholesky decomposition of matrix
function getChol(mat::Array)
    #get spd here, then just do chol(R)'
    #U,E,V=svd(mat)
    return Matrix(chol(nearestSPD(mat)))
end

# Calculation for updating c
function updateC(m::AbstractModel,c_in::Float64, acpt_in::Float64, trgt_in::Float64)
    c_out=c_in*(0.95+0.1*exp(16*(acpt_in-trgt_in))/(1+exp(16*(acpt_in-trgt_in))))
    m<=Setting(:tpf_c,c_out)
    return c_out
end

################
function correctAndResample(φ_new::Float64, φ_old::Float64, yt::Array, perror::Array,density_arr::Array,weights::Array, s_up::Array, ε_up::Array, H::Array, numParticles::Int64; initialize=0)

    # Calculate initial weights
    for n=1:numParticles
        density_arr[n]=density(φ_new, φ_old, yt, perror[:,n], H, initialize=initialize)
    end
   
    # Normalize weights
    weights = (density_arr.*weights)./mean(density_arr.*weights)
    # Resampling
    id = multinomial_resampling(weights)
    #rintln("IDs: $id")
    s_up = s_up[:,id]
    ε_up = ε_up[:,id]
    
    # Reset weights to ones
    weights = ones(numParticles)

    # Calculate likelihood
    loglik = log(mean(density_arr.*weights))

    return loglik, weights, s_up, ε_up
end
##############3
