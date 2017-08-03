using Roots
function tpf(m::AbstractModel, yy::Array, system::System{Float64}, s0::Array{Float64}, P0::Array; verbose::Symbol=:low, include_presample::Bool=true)
    # s0 is 8xn_particles
    # P0 is solution to discrete lyapunov equation for Φ and R*S2*R'
    # yy is data matrix

    #--------------------------------------------------------------
    # Set Parameters of Algorithm
    #--------------------------------------------------------------

    # Store model parameters
    RRR    = system.transition.RRR
    TTT    = system.transition.TTT
    EE     = system.measurement.EE
    DD     = system.measurement.DD
    ZZ     = system.measurement.ZZ
    QQ     = system.measurement.QQ    
    sqrtS2 = RRR*get_chol(QQ)'

    # Get tuning parameters from the model
    rstar              = get_setting(m, :tpf_rstar)
    c                  = get_setting(m, :tpf_c)
    acpt_rate          = get_setting(m, :tpf_acpt_rate)
    target             = get_setting(m, :tpf_target)
    N_MH               = get_setting(m, :tpf_n_mh_simulations)
    n_particles        = get_setting(m, :tpf_n_particles)
    deterministic      = get_setting(m, :tpf_deterministic)
    xtol               = get_setting(m, :tpf_x_tolerance) # Tolerance of fzero
    parallel           = get_setting(m, :use_parallel_workers) # Get setting of parallelization
    mutation_rand_mat  = get_setting(m, :tpf_rand_mat)
    
    # Determine presampling periods
    n_presample_periods = (include_presample) ? 0 : get_setting(m, :n_presample_periods)
    
    # End time (last period)
    T = size(yy,2)
    # Size of covariance matrix
    n_errors = size(QQ,1)
    # Number of states
    n_states = size(ZZ,2)
    
    # Initialization
    lik         = zeros(T)
    Neff        = zeros(T)
    len_phis    = ones(T)
    times       = zeros(T)
    weights     = ones(n_particles)
    density_arr = zeros(n_particles)
    cov_s       = eye(n_errors)
    
    #resampling_ids = zeros(3*T,n_particles)
    #ids_i = 1

    if deterministic
        # Testing: Generate pre-defined random shocks for s and ε
        s_lag_tempered_rand_mat = randn(n_states,n_particles)
        ε_rand_mat = randn(n_errors, n_particles)
        path = dirname(@__FILE__)  
        h5open("$path/../../test/reference/matricesForTPF.h5","w") do file
            write(file,"s_lag_tempered_rand_mat",s_lag_tempered_rand_mat)
            write(file,"eps_rand_mat",ε_rand_mat)
        end
    else
        # Draw initial particles from the distribution of s₀: N(s₀, P₀) 
        s_lag_tempered_rand_mat = randn(n_states,n_particles)
    end

    ###change back to get_chol later!!
	s_lag_tempered = repmat(s0,1,n_particles) + Matrix(chol(P0))'*s_lag_tempered_rand_mat
    times = zeros(T)

    for t=1:T

        tic()
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            println("============================================================")
            @show t
        end
        
        #--------------------------------------------------------------
        # Initialize Algorithm: First Tempering Step
        #--------------------------------------------------------------
        y_t = yy[:,t]
        
        # Remove rows/columns of series with NaN values
        nonmissing = !isnan(y_t)
        y_t         = y_t[nonmissing]        
        ZZ_t       = ZZ[nonmissing,:]
        DD_t       = DD[nonmissing]
        EE_t       = EE[nonmissing,nonmissing]
        QQ_t       = QQ[nonmissing,nonmissing]
        RRR_t      = RRR[:,nonmissing]
        sqrtS2_t   = RRR_t*get_chol(QQ_t)'
        n_errors_t = length(y_t)
        ε          = zeros(n_errors_t)
            
        if !deterministic # When not testing, want a new random ε each time 
            ε_rand_mat = randn(n_errors_t, n_particles)
        end

        # TEMPORARY
        ε_rand_mat = randn(n_errors_t, n_particles)

        # Forecast forward one time step
        s_t_nontempered = TTT*s_lag_tempered + sqrtS2_t*ε_rand_mat
        
        # Error for each particle
        p_error = repmat(y_t - DD_t, 1, n_particles) - ZZ_t*s_t_nontempered

        # Solve for initial tempering parameter ϕ_1
        if !deterministic
            init_Ineff_func(φ) = ineff_func(φ, 2.0*pi, y_t, p_error, EE_t, initialize=1)-rstar
            φ_1 = fzero(init_Ineff_func,0.00000001, 1.0, xtol=xtol)
        else 
            φ_1 = 0.25
        end
        
         if VERBOSITY[verbose] >= VERBOSITY[:low]
            @show φ_1 
            println("------------------------------")
        end
              
        # Update weights array and resample particles
        loglik, weights, s_lag_tempered, ε, id = correct_and_resample(φ_1,0.0,y_t,p_error,density_arr,weights,s_lag_tempered,ε_rand_mat,EE_t,n_particles,deterministic,initialize=1)
        #resampling_ids[ids_i,:] = id
        #ids_i += 1

        # Update likelihood
        lik[t] += loglik
        
        # Tempering Initialization
        count = 2 # Accounts for initialization and final mutation
        φ_old = φ_1

        # First propagation
        s_t_nontempered = TTT*s_lag_tempered + sqrtS2_t*ε
        p_error = repmat(y_t - DD_t, 1, n_particles) - ZZ_t*s_t_nontempered         
        
        if !deterministic
            println("You're not deterministic!")
            check_ineff = ineff_func(1.0, φ_1, y_t, p_error, EE_t)         
        else
            check_ineff = rstar + 1
        end

        if VERBOSITY[verbose] >= VERBOSITY[:high]
            @show check_ineff
        end

        #--------------------------------------------------------------
        # Main Algorithm
        #--------------------------------------------------------------
        while check_ineff > rstar

            # Define inefficiency function
            init_ineff_func(φ) = ineff_func(φ, φ_old, y_t, p_error, EE_t) - rstar
            φ_interval = [φ_old, 1.0]
            fphi_interval = [init_ineff_func(φ_old) init_ineff_func(1.0)]

            count += 1

            # Check solution exists within interval
            if prod(sign(fphi_interval)) == -1 || deterministic
                
                if deterministic
                    φ_new = 0.5
                else
                    # Set φ_new to the solution of the inefficiency function over interval
                    φ_new = fzero(init_ineff_func, φ_interval, xtol=xtol)
                    check_ineff = ineff_func(1.0, φ_new, y_t, p_error, EE_t, initialize=0)
                end
               
                if VERBOSITY[verbose] >= VERBOSITY[:low]
                    @show φ_new
                end

                # Update weights array and resample particles
                loglik, weights, s_lag_tempered, ε, id = correct_and_resample(φ_new, φ_old, y_t, p_error, density_arr, weights, s_lag_tempered, ε, EE_t, n_particles, deterministic, initialize=0)
                #resampling_ids[ids_i,:] = id
                #ids_i += 1

                # Update likelihood
                lik[t] += loglik
                
                # Update value for c
                c = update_c!(m, c, acpt_rate, target)
                
                if VERBOSITY[verbose] >= VERBOSITY[:low]
                    @show c
                    println("------------------------------")
                end
                                
                # Update covariance matrix
                μ = mean(ε,2)
                if !deterministic
                    cov_s = (1/n_particles)*(ε - repmat(μ,1,n_particles))*(ε - repmat(μ,1,n_particles))'
                    if !isposdef(cov_s) cov_s = diagm(diag(cov_s)) end
                end
                
                # Mutation Step
                acpt_vec = zeros(n_particles)
                print("Mutation ")        
                tic()

                if parallel
                    print("(in parallel) ")
                    c = get_setting(m,:tpf_c)
                    N_MH = get_setting(m,:tpf_n_mh_simulations)
                    deterministic = get_setting(m,:tpf_deterministic)
                    #out = pmap(i->mutation(system,y_t,s_lag_tempered[:,i],ε[:,i],c, N_MH,deterministic,cov_s,nonmissing,mutation_rand_mat), 1:n_particles)
                    out = @sync @parallel (hcat) for i=1:n_particles
                        mutation(system,y_t,s_lag_tempered[:,i],ε[:,i],c,N_MH,deterministic,cov_s,nonmissing,mutation_rand_mat)
                    end
                else 
                    print("(not parallel) ")
                    out = [mutation(system,y_t,s_lag_tempered[:,i],ε[:,i],c,N_MH,deterministic,cov_s,nonmissing,mutation_rand_mat) for i=1:n_particles]
                end
                times[t] = toc()                

                for i = 1:n_particles
                    s_t_nontempered[:,i] = out[i][1]
                    ε[:,i] = out[i][2]
                    acpt_vec[i] = out[i][3]
                end

                # Calculate average acceptance rate
                acpt_rate = mean(acpt_vec)

                # Get error for all particles
                p_error = repmat(y_t-DD_t, 1, n_particles)-ZZ_t*s_t_nontempered
                φ_old = φ_new
                len_phis[t] += 1

            # If no solution exists within interval, set inefficiency to rstar
            else 
                if VERBOSITY[verbose] >= VERBOSITY[:high]
                    println("No solution in interval.")
                end
                check_ineff = rstar
            end
            gc()

            # With phi schedule, leave while loop after one iteration, thus set check_ineff=0
            if deterministic
                check_ineff = 0
            end

            if VERBOSITY[verbose] >= VERBOSITY[:high]
                @show check_ineff
            end
        end

        if VERBOSITY[verbose] >= VERBOSITY[:high]
            println("Out of main while-loop.")
        end
        
        #--------------------------------------------------------------
        # Last Stage of Algorithm: φ_new=1
        #--------------------------------------------------------------
        φ_new = 1.0

        # Update weights array and resample particles.
        loglik, weights, s_lag_tempered, ε, id = correct_and_resample(φ_new,φ_old,y_t,p_error,density_arr,weights,s_lag_tempered,ε,EE_t,n_particles,deterministic,initialize=0)
        #resampling_ids[ids_i,:] = id
        #ids_i += 1

        # Update likelihood
        lik[t] += loglik

        c = update_c!(m, c, acpt_rate, target)
        μ = mean(ε,2)

        if !deterministic
            cov_s = (1/n_particles)*(ε - repmat(μ,1,n_particles))*(ε - repmat(μ,1,n_particles))'
            if !isposdef(cov_s) cov_s = diagm(diag(cov_s)) end
        end
        
        # Final round of mutation
        acpt_vec = zeros(n_particles)

        if parallel
            c = get_setting(m,:tpf_c)
            N_MH = get_setting(m,:tpf_n_mh_simulations)
            deterministic = get_setting(m,:tpf_deterministic)

            # out = pmap(i -> mutation(system,y_t,s_lag_tempered[:,i],ε[:,i],c,N_MH,deterministic,cov_s,nonmissing,mutation_rand_mat), 1:n_particles)
            out = @sync @parallel (hcat) for i=1:n_particles
                mutation(system,y_t,s_lag_tempered[:,i],ε[:,i],c,N_MH,deterministic,cov_s,nonmissing,mutation_rand_mat)
            end
        else 
            out = [mutation(system,y_t,s_lag_tempered[:,i],ε[:,i],c,N_MH,deterministic,cov_s,nonmissing,mutation_rand_mat) for i=1:n_particles]
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
        print("Completion of one period ")
        gc()
        toc()
    end

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("=============================================")
    end
#= 
    if deterministic
        h5open("$path/../../test/reference/resampled_ids.h5","w") do f
            write(f, "resampling_ids",resampling_ids)
        end
    end
=#
    # Return vector of likelihood indexed by time step and Neff
    return Neff[n_presample_periods + 1:end], lik[n_presample_periods + 1:end], times
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
update_c!(m::AbstractModel, c_in::Float64, acpt_in::Float64, target_in::Float64)
```
Update value of c by expression that is function of the target and mean acceptance rates.
Returns the new c, in addition to storing it in the model settings.

"""
function update_c!(m::AbstractModel,c_in::Float64, acpt_in::Float64, target_in::Float64)
    c_out = c_in*(0.95 + 0.1*exp(16*(acpt_in - target_in))/(1 + exp(16*(acpt_in - target_in))))
    m <= Setting(:tpf_c, c_out)
    return c_out
end

"""
```
correct_and_resample(φ_new::Float64, φ_old::Float64, y_t::Array, p_error::Array,density_arr::Array,weights::Array, s_lag_tempered::Array, ε::Array, EE::Array, n_particles::Int64; initialize::Int64=0)
```
Calculate densities, normalize and reset weights, call multinomial resampling, update state and error vectors,reset error vectors to 1,and calculate new log likelihood.
Returns log likelihood, weight, state, and ε vectors.

"""
function correct_and_resample(φ_new::Float64, φ_old::Float64, y_t::Array{Float64,1}, p_error::Array{Float64,2}, density_arr::Array{Float64,1}, weights::Array{Float64,1}, s_lag_tempered::Array{Float64,2}, ε::Array{Float64,2}, EE::Array{Float64,2}, n_particles::Int64, deterministic::Bool; initialize::Int64=0)
    # Calculate initial weights
    for n=1:n_particles
        density_arr[n]=density(φ_new, φ_old, y_t, p_error[:,n], EE, initialize=initialize)
    end   

    # Normalize weights
    ##SHOULD THIS BE SUM INSTEAD OF MEAN?
    weights = (density_arr.*weights)./mean(density_arr.*weights)
    
    # Resampling
    if !deterministic
        id = multinomial_resampling(weights)
    else
        id = seeded_multinomial_resampling(weights)
    end
    
    s_lag_tempered = s_lag_tempered[:,id]
    ε = ε[:,id]

    # Reset weights to ones
    weights = ones(n_particles)

    # Calculate likelihood
    loglik = log(mean(density_arr.*weights))
    
    return loglik, weights, s_lag_tempered, ε, id
end

function zlb_regime_indices{S<:AbstractFloat}(m::AbstractModel{S},data::Matrix{S})
    # Make sure the data matrix has all time periods when passing in or this won't work
    T = size(data,2)
    if n_anticipated_shocks(m) > 0
        regime_inds = Vec{Range{Int64}}(2)
        regime_inds[1] = 1 : index_zlb_start(m) - 1
        regime_inds[2] = index_zlb_start(m) : T 
    else 
        regime_inds = Range{Int64}[1:T]
    end
end

function zlb_regime_matrices{S<:AbstractFloat}(m::AbstractModel{S},system::System{S})
    if !all(x -> x==0, system[:MM])
        error("Previously this error said Kalman filter and smoothers not implemented for nonzero MM however i'm not sure if this still applies to the TPF")
    end
    
    if n_anticipated_shocks(m) > 0
        n_regimes = 2
        
        shock_inds = inds_shocks_no_ant(m)
        QQ_ZLB = zeros(size(QQ_ZLB))
        QQ_preZLB[shock_inds, shock_inds] = QQ_ZLB[shock_inds,shock_inds]
        QQs = Matrix{S}[QQ_preZLB,QQ_ZLB]
    else 
        n_regimes = 1
        QQs = Matrix{S}[system[:QQ]]
    end
    TTTs = fill(system[:TTT], n_regimes)
    RRRs = fill(system[:RRR], n_regimes)
    CCCs = fill(system[:CCC], n_regimes)
    ZZs  = fill(system[:ZZ],  n_regimes)
    DDs  = fill(system[:DD],  n_regimes)
    EEs  = fill(system[:EE],  n_regimes)

    return TTTs, RRRs, CCCs, QQs, ZZs, DDs, EEs
end


nothing
