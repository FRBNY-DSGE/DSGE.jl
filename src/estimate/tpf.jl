using DSGE
using Roots
function tpf(m::AbstractModel, yy::Array, s0::Array{Float64}, P0::Array, A, B, H, R, S2, Φ)
    # s0 is 8xnum_particles
    # P0
    # yy is data matrix

#--------------------------------------------------------------
# Set Parameters of Algorithm
#--------------------------------------------------------------

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
    sqrtS2=R*get_chol(S2)'

    path = dirname(@__FILE__)


    # Get tuning parameters from the model
  
##### WILL ALSO NEED TO BRING THIS BACK LATER
    #= rstar = get_setting(m,:tpf_rstar)
    c = get_setting(m,:tpf_c)
    acpt_rate = get_setting(m,:tpf_acpt_rate)
    trgt = get_setting(m,:tpf_trgt)
    N_MH = get_setting(m,:tpf_N_MH)
    num_particles = get_setting(m,:tpf_num_particles)
    =#
#### BELOW IS HARD CODED, TO BE REPLACED BY ABOVE
    rstar= 2.0
    acpt_rate= 0.5
    trgt= 0.25
    c = 0.1
    N_MH= 2
    num_particles=100

    # End time (last period)
    T = size(yy,2)
    # Size of covariance matrix
    num_errors = size(S2,1)
    #number of states
    num_states = size(B, 2)
    
    # Likelihood matrix
    lik = zeros(T)
    Neff=zeros(T)
    len_phis = ones(T)
    weights = ones(num_particles)
    density_arr=zeros(num_particles)
    ε_up = zeros(num_errors)

    # Initialize first state
    # STORING FOR TESTING - RANDOM
    s_up_rand_mat = randn(num_states,num_particles)
    ε_rand_mat = randn(num_errors, num_particles)
    h5open("$path/../../test/reference/matricesForTPF.h5","w") do file
        write(file,"s_up_rand_mat",s_up_rand_mat)
        write(file,"ε_rand_mat",ε_rand_mat)
    end
    ###
    s_up = repmat(s0,1,num_particles) + get_chol(P0)'*s_up_rand_mat
    
    for t=1:T

        #--------------------------------------------------------------
        # Initialize Algorithm: First Tempering Step
        #--------------------------------------------------------------

        yt = yy[:,t]
        # RANDOM
        #ε = randn(num_errors, num_particles)
        s_fore = Φ*s_up + sqrtS2*ε_rand_mat

        # Error for each particle
        perror = repmat(yt-A,1,num_particles)-B*s_fore

        # Solve for initial tempering parameter ϕ_1
        init_Ineff_func(φ) = ineff_func(φ, 2.0*pi, yt, perror, H, initialize=1)-rstar
        φ_1 = fzero(init_Ineff_func,0.000001, 1.0)

        # Update weights array and resample particles
        loglik, weights, s_up, ε_up = correct_and_resample(φ_1,0.0,yt,perror,density_arr,weights,s_up,ε,H,num_particles,initialize=1)
        # Update likelihood
        lik[t]=lik[t]+loglik

        # Tempering Initialization
        count = 2 # Accounts for initialization and final mutation
        φ_old = φ_1

        # First propagation
        s_fore = Φ*s_up + sqrtS2*ε_up
        perror = repmat(yt-A, 1, num_particles) - B*s_fore        
        check_ineff=ineff_func(1.0, φ_1, yt, perror, H)         

        #--------------------------------------------------------------
        # Main Algorithm
        #--------------------------------------------------------------
        while (check_ineff>rstar)

            count += 1

            # Define inefficiency function
            init_ineff_func(φ) = ineff_func(φ, φ_old, yt, perror,H)-rstar
            φ_interval = [φ_old, 1.0]
            fphi_interval = [init_ineff_func(φ_old) init_ineff_func(1.0)]
            
            # Check solution exists within interval
            if prod(sign(fphi_interval))==-1
               
                # Set φ_new to the solution of the inefficiency function over interval
                φ_new = fzero(init_ineff_func, φ_interval)
                check_ineff = ineff_func(1.0, φ_new, yt, perror, H, initialize=0)
               
                # Update weights array and resample particles
                loglik, weights, s_up, ε_up = correct_and_resample(φ_new,φ_old,yt,perror,density_arr,weights,s_up,ε_up,H,num_particles,initialize=0)

                # Update likelihood
                lik[t]=lik[t]+loglik

                # Update value for c
                c = update_c(m,c,acpt_rate,trgt)
                
                # Update covariance matrix
                μ = mean(ε_up,2)
                cov_s = (1/num_particles)*(ε_up-repmat(μ,1,num_particles))*(ε_up - repmat(μ,1,num_particles))'
 
                # Mutation Step
                acpt_vec=zeros(num_particles)
                for i = 1:num_particles
                    ind_s, ind_ε, ind_acpt = mutation(m,yt,s_up[:,i],ε_up[:,i],A,B,cov_s,Φ,H,sqrtS2,S2,N_MH)
                    s_fore[:,i] = ind_s
                    ε_up[:,i] = ind_ε
                    acpt_vec[i] = ind_acpt
                end

                # Calculate average acceptance rate
                acpt_rate = mean(acpt_vec)

                # Get error for all particles
                perror = repmat(yt-A, 1,num_particles)-B*s_fore
                φ_old = φ_new
                len_phis[t]+=1

            # If no solution exists within interval, set inefficiency to rstar
            else 
                check_ineff = rstar
            end
        end

        #--------------------------------------------------------------
        # Last Stage of Algorithm: φ_new=1
        #--------------------------------------------------------------
        φ_new = 1.0

        # Update weights array and resample particles.
        loglik, weights, s_up, ε_up = correct_and_resample(φ_new,φ_old,yt,perror,density_arr,weights,s_up,ε_up,H,num_particles,initialize=0)
        # Update likelihood
        lik[t]=lik[t]+loglik
        
        c = update_c(m,c,acpt_rate,trgt)

        μ = mean(ε_up,2)
        cov_s = (1/num_particles)*(ε_up-repmat(μ,1,num_particles))*(ε_up - repmat(μ,1,num_particles))'
       
        # Final round of mutation
        acpt_vec=zeros(num_particles)
        for i=1:num_particles
            ind_s, ind_ε, ind_acpt = mutation(m, yt, s_up[:,i], ε_up[:,i], A, B, cov_s, Φ, H, sqrtS2,S2,N_MH)
            s_fore[:,i] = ind_s
            ε_up[:,i] = ind_ε
            acpt_vec[i] = ind_acpt 
        end

        # Store for next time iteration
        acpt_rate = mean(acpt_vec)
        Neff[t] = (num_particles^2)/sum(weights.^2)
        s_up = s_fore
    end
    # Return vector of likelihood indexed by time step and Neff
    return Neff, lik
end


"""
```
get_chol(mat::Aray)
```
Calculate and return the Cholesky of a matrix.

"""
function get_chol(mat::Array)
    return Matrix(chol(nearestSPD(mat)))
end

"""
```
update_c(m::AbstractModel, c_in::Float64, acpt_in::Float64, trgt_in::Float64)
```
Update value of c by expression that is function of the target and mean acceptance rates.
Returns the new c, in addition to storing it in the model settings.

"""
function update_c(m::AbstractModel,c_in::Float64, acpt_in::Float64, trgt_in::Float64)
    c_out=c_in*(0.95+0.1*exp(16*(acpt_in-trgt_in))/(1+exp(16*(acpt_in-trgt_in))))
    m<=Setting(:tpf_c,c_out)
    return c_out
end

"""
```
correct_and_resample(φ_new::Float64, φ_old::Float64, yt::Array, perror::Array,density_arr::Array,weights::Array, s_up::Array, ε_up::Array, H::Array, num_particles::Int64; initialize::Int64=0)
```
Calculate densities, normalize and reset weights, call multinomial resampling, update state and error vectors,reset error vectors to 1,and calculate new log likelihood.
Returns log likelihood, weight, state, and ε vectors.

"""
function correct_and_resample(φ_new::Float64, φ_old::Float64, yt::Array, perror::Array, density_arr::Array, weights::Array, s_up::Array, ε_up::Array, H::Array, num_particles::Int64; initialize::Int64=0)

    # Calculate initial weights
    for n=1:num_particles
        density_arr[n]=density(φ_new, φ_old, yt, perror[:,n], H, initialize=initialize)
    end   

    # Normalize weights
    weights = (density_arr.*weights)./mean(density_arr.*weights)

    # Resampling
    id = multinomial_resampling(weights)
    s_up = s_up[:,id]
    ε_up = ε_up[:,id]
    
    # Reset weights to ones
    weights = ones(num_particles)

    # Calculate likelihood
    loglik = log(mean(density_arr.*weights))

    return loglik, weights, s_up, ε_up
end

nothing
