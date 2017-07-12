using DSGE
using Roots
function tpf(m::AbstractModel, yy::Array, s0::Array{Float64}, P0::Array; testing::Bool=false,verbose::Symbol=:low)
    # s0 is 8xn_particles
    # P0
    # yy is data matrix

    #--------------------------------------------------------------
    # Set Parameters of Algorithm
    #--------------------------------------------------------------

    # Compute system and store model parameters
    if !testing
        # Non-testing mode. Use actual matrices computed from model system.
        sys=compute_system(m)
        R=sys.transition.RRR
        Φ=sys.transition.TTT
        H=sys.measurement.EE
        A=sys.measurement.DD
        B=sys.measurement.ZZ
        S2=sys.measurement.QQ
 
        m<=Setting(:DD,A)
        m<=Setting(:ZZ,B)
        m<=Setting(:RRR,R)
        m<=Setting(:TTT,Φ)
        m<=Setting(:EE,H)
        m<=Setting(:tpf_S2,S2)
    else   
        # Testing Mode. Read in matrices from Schorfheide Matlab code.
        A = get_setting(m,:DD)
        B = get_setting(m,:ZZ)
        R = get_setting(m,:RRR)
        Φ = get_setting(m,:TTT)
        H = get_setting(m,:EE)
        S2 = get_setting(m,:tpf_S2)
    end   
    
    sqrtS2=R*get_chol(S2)'
    m<=Setting(:tpf_sqrtS2,sqrtS2)
    
    # Get tuning parameters from the model
    rstar = get_setting(m,:tpf_rstar)
    c = get_setting(m,:tpf_c)
    acpt_rate = get_setting(m,:tpf_acpt_rate)
    trgt = get_setting(m,:tpf_trgt)
    N_MH = get_setting(m,:tpf_N_MH)
    n_particles = get_setting(m,:tpf_n_particles)
    xtol = get_setting(m, :x_tolerance)
    parallel = get_setting(m,:use_parallel_workers)
    path = dirname(@__FILE__)  
    
    # Testing mode: Generate pre-defined 3x1 vector of random movements for mutation step. 
    rand_mat = randn(3,1)
    
    if testing
        h5open("$path/../../test/reference/mutationRandomMatrix.h5","w") do file
            write(file,"rand_mat",rand_mat)
        end
    end

    # End time (last period)
    T = size(yy,2)
    # Size of covariance matrix
    n_errors = size(S2,1)
    # Number of states
    n_states = size(B,2)
    
    # Initialization
    lik = zeros(T)
    Neff=zeros(T)
    len_phis = ones(T)
    weights = ones(n_particles)
    density_arr=zeros(n_particles)
    ε = zeros(n_errors)

    # Initialize first state

    ### Testing: Generated pre-defined random shocks for s and ε
    if testing
        s_lag_tempered_rand_mat = randn(n_states,n_particles)
        ε_rand_mat = randn(n_errors, n_particles)
        h5open("$path/../../test/reference/matricesForTPF$n_particles.h5","w") do file
            write(file,"s_lag_tempered_rand_mat",s_lag_tempered_rand_mat)
            write(file,"eps_rand_mat",ε_rand_mat)
        end
    end
    ###

    #Draw initial particles from the distribution of s₀: N(s₀, P₀) 
    if !testing #if not testing, draw a new random term. 
        s_lag_tempered_rand_mat = randn(n_states,n_particles)
    end

    s_lag_tempered = repmat(s0,1,n_particles) + get_chol(P0)'*s_lag_tempered_rand_mat

    for t=1:T
        tic()
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            @show t
        end
        
        #--------------------------------------------------------------
        # Initialize Algorithm: First Tempering Step
        #--------------------------------------------------------------
        yt = yy[:,t]
        
        if !testing # When not testing, want a new random ε each time 
            ε_rand_mat = randn(n_errors, n_particles)
        end

        # Forecast forward one time step
        s_t_nontempered = Φ*s_lag_tempered + sqrtS2*ε_rand_mat
        
        # Error for each particle
        perror = repmat(yt-A,1,n_particles)-B*s_t_nontempered
        
        # Solve for initial tempering parameter ϕ_1
        if !testing
            init_Ineff_func(φ) = ineff_func(φ, 2.0*pi, yt, perror, H, initialize=1)-rstar
            φ_1 = fzero(init_Ineff_func,0.000001, 1.0, xtol=xtol)
        else 
            φ_1 = 0.25
        end
        
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            @show φ_1 
        end
              
        # Update weights array and resample particles
        # While it might seem weird that we're updating s_lag_tempered and not s_t_nontempered, we are actually just resampling and we only need the lagged states for future calculations.
        loglik, weights, s_lag_tempered, ε = correct_and_resample(φ_1,0.0,yt,perror,density_arr,weights,s_lag_tempered,ε_rand_mat,H,n_particles,testing,initialize=1)
        
        # Update likelihood
        lik[t]=lik[t]+loglik

        # Tempering Initialization
        count = 2 # Accounts for initialization and final mutation
        φ_old = φ_1

        # First propagation
        
        #turns out these lines are actually important! figure out later why
        s_t_nontempered = Φ*s_lag_tempered + sqrtS2*ε
        perror = repmat(yt-A, 1, n_particles) - B*s_t_nontempered         
        
        if !testing
            @show "You're not testing!"
            check_ineff=ineff_func(1.0, φ_1, yt, perror, H)         
        else
            check_ineff=rstar+1
        end

        if VERBOSITY[verbose] >= VERBOSITY[:low]
            @show check_ineff
        end
        #--------------------------------------------------------------
        # Main Algorithm
        #--------------------------------------------------------------
        while (check_ineff>rstar)
            # Define inefficiency function
            init_ineff_func(φ) = ineff_func(φ, φ_old, yt, perror,H)-rstar
            φ_interval = [φ_old, 1.0]
            fphi_interval = [init_ineff_func(φ_old) init_ineff_func(1.0)]
            count += 1

            # Check solution exists within interval
            if prod(sign(fphi_interval))== -1 || testing
                
                if !testing
                    # Set φ_new to the solution of the inefficiency function over interval
                    println("Time of fphi inside while loop")
                    tic()
                    φ_new = fzero(init_ineff_func, φ_interval, xtol=xtol)
                    toc()
                    check_ineff = ineff_func(1.0, φ_new, yt, perror, H, initialize=0)
                else
                    φ_new = 0.5
                end
               
                if VERBOSITY[verbose] >= VERBOSITY[:low]
                    @show φ_new
                end

                # Update weights array and resample particles
                loglik, weights, s_lag_tempered, ε = correct_and_resample(φ_new,φ_old,yt,perror,density_arr,weights,s_lag_tempered,ε,H,n_particles, testing, initialize=0)

                # Update likelihood
                lik[t] = lik[t] + loglik
                                
                # Update value for c
                c = update_c(m,c,acpt_rate,trgt)
                
                if VERBOSITY[verbose] >= VERBOSITY[:low]
                    @show c
                end
                                
                # Update covariance matrix
                μ = mean(ε,2)
                cov_s = (1/n_particles)*(ε-repmat(μ,1,n_particles))*(ε - repmat(μ,1,n_particles))'
                
                if !isposdef(cov_s)
                    cov_s = diagm(diag(cov_s))
                end
                
                # Mutation Step
                acpt_vec=zeros(n_particles)
                @show "Mutation time..."
                tic()
                if parallel
                    @show "In parallel:"
                    out = pmap(i->mutation(m,yt,s_lag_tempered[:,i],ε[:,i],cov_s,rand_mat,testing), 1:n_particles)
                else 
                    @show "Not in parallel:"
                    out = [mutation(m,yt,s_lag_tempered[:,i],ε[:,i],cov_s,rand_mat,testing) for i=1:n_particles]
                end
                toc()                
                
                for i = 1:n_particles
                    s_t_nontempered[:,i] = out[i][1]
                    ε[:,i] = out[i][2]
                    acpt_vec[i] = out[i][3]
                end

                # Calculate average acceptance rate
                acpt_rate = mean(acpt_vec)

                # Get error for all particles
                perror = repmat(yt-A, 1,n_particles)-B*s_t_nontempered
                φ_old = φ_new
                len_phis[t]+=1

            # If no solution exists within interval, set inefficiency to rstar
            else 
                if VERBOSITY[verbose] >= VERBOSITY[:high]
                    println("No solution in interval.")
                end
                check_ineff = rstar
            end
            #For testing with phi schedule, leave while loop after one iteration, thus set check_ineff=0
            if testing
               check_ineff=0
            end

            if VERBOSITY[verbose] >= VERBOSITY[:low]
                @show check_ineff
            end
        end

        if VERBOSITY[verbose] >= VERBOSITY[:low]
            println("Out of main while-loop.")
        end
        
        #--------------------------------------------------------------
        # Last Stage of Algorithm: φ_new=1
        #--------------------------------------------------------------
        φ_new = 1.0

        # Update weights array and resample particles.
        loglik, weights, s_lag_tempered, ε = correct_and_resample(φ_new,φ_old,yt,perror,density_arr,weights,s_lag_tempered,ε,H,n_particles,testing,initialize=0)

       # Update likelihood
        lik[t]=lik[t]+loglik

        c = update_c(m,c,acpt_rate,trgt)
        μ = mean(ε,2)

        cov_s = (1/n_particles)*(ε-repmat(μ,1,n_particles))*(ε - repmat(μ,1,n_particles))'
                
        if isposdef(cov_s)== false 
            cov_s = diagm(diag(cov_s))
        end
        
        # Final round of mutation
        acpt_vec=zeros(n_particles)
  
        if parallel
            out = pmap(i -> mutation(m,yt,s_lag_tempered[:,i],ε[:,i],cov_s,rand_mat,testing), 1:n_particles)
        else 
            out = [mutation(m, yt, s_lag_tempered[:,i], ε[:,i], cov_s, rand_mat, testing) for i=1:n_particles]
        end
                
        for i = 1:n_particles
            s_t_nontempered[:,i] = out[i][1]
            ε[:,i] = out[i][2]
            acpt_vec[i] = out[i][3]
        end
        
        # Store for next time iteration
        acpt_rate = mean(acpt_vec)

        Neff[t] = (n_particles^2)/sum(weights.^2)
        s_lag_tempered = s_t_nontempered
        println("Time of entire time loop")
        toc()
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
function get_chol(mat::Array{Float64,2})
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
correct_and_resample(φ_new::Float64, φ_old::Float64, yt::Array, perror::Array,density_arr::Array,weights::Array, s_lag_tempered::Array, ε::Array, H::Array, n_particles::Int64; initialize::Int64=0)
```
Calculate densities, normalize and reset weights, call multinomial resampling, update state and error vectors,reset error vectors to 1,and calculate new log likelihood.
Returns log likelihood, weight, state, and ε vectors.

"""
function correct_and_resample(φ_new::Float64, φ_old::Float64, yt::Array{Float64,1}, perror::Array{Float64,2}, density_arr::Array{Float64,1}, weights::Array{Float64,1}, s_lag_tempered::Array{Float64,2}, ε::Array{Float64,2}, H::Array{Float64,2}, n_particles::Int64, testing::Bool; initialize::Int64=0)

    # Calculate initial weights
    path = dirname(@__FILE__)
    for n=1:n_particles
        density_arr[n]=density(φ_new, φ_old, yt, perror[:,n], H, initialize=initialize)
    end   
        
    # Normalize weights
    weights = (density_arr.*weights)./mean(density_arr.*weights)

    # Resampling
    if !testing
        id = multinomial_resampling(weights)
    else
        id = seeded_multinomial_resampling(weights)
    end
    # @show id
    # fid1 = open("$path/../../test/reference/random_ids.csv","a")
    # for i in id
    #     print(fid1,i,",")
    # end    
    # print(fid1,"\n")
    # close(fid1)

    s_lag_tempered = s_lag_tempered[:,id]
    ε = ε[:,id]

    # Reset weights to ones
    weights = ones(n_particles)

    # Calculate likelihood
    loglik = log(mean(density_arr.*weights))
    
    return loglik, weights, s_lag_tempered, ε
end

nothing
