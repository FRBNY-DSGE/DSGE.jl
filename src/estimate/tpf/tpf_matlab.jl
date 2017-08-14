using DSGE
function tpf_matlab{S<:AbstractFloat}(m::AbstractModel, yy::Array, system::System{S},
    s0::Array{S}, P0::Array; verbose::Symbol=:low, include_presample::Bool=true)
    # s0 is 8xn_particles
    # P0 is solution to discrete lyapunov equation for Φ and R*S2*R'
    # yy is data matrix

    #--------------------------------------------------------------
    # Set Parameters of Algorithm
    #--------------------------------------------------------------

    # Unpack system
    RRR = system[:RRR]
    TTT = system[:TTT]
    EE  = system[:EE]
    DD  = system[:DD]
    ZZ  = system[:ZZ]
    QQ  = system[:QQ]    
    
    # Get tuning parameters from the model
    rstar              = get_setting(m, :tpf_rstar)
    c                  = get_setting(m, :tpf_c)
    accept_rate        = get_setting(m, :tpf_accept_rate)
    target             = get_setting(m, :tpf_target)
    N_MH               = get_setting(m, :tpf_n_mh_simulations)
    n_particles        = get_setting(m, :tpf_n_particles)
    deterministic      = get_setting(m, :tpf_deterministic)
    xtol               = get_setting(m, :tpf_x_tolerance)
    parallel           = get_setting(m, :use_parallel_workers)
    mutation_rand_mat  = get_setting(m, :tpf_rand_mat)    

    # Determine presampling periods
    n_presample_periods = (include_presample) ? 0 : get_setting(m, :n_presample_periods)
    
    # End time (last period)
    T = size(yy,2) 
    # Size of covariance matrix
    n_observables = size(QQ,1)
    # Number of states
    n_states = size(ZZ,2)
    
    # Initialization
    lik                 = zeros(T)
    Neff                = zeros(T)
    len_phis            = ones(T)
    times               = zeros(T)
    weights             = ones(n_particles)
    incremental_weights = zeros(n_particles)
    times               = zeros(T)
    
    # resampling_ids = zeros(3*T,n_particles)
    # ids_i = 1
    id_cur = 1
    id_m = 1
    s_id = 1 

    # TEMPORARY: 
    matlab_resampled_ids = h5read("/data/dsge_data_dir/dsgejl/interns2017/matlab_resampled_ids.h5",
                                  "resampled_ids")
    matlab_resampled_ids = round(Int64, matlab_resampled_ids)
    mutation_rand_mats1 = h5read("/data/dsge_data_dir/dsgejl/interns2017/mutation_rand_mats.h5",
                                 "mutation_rand_mats1")
    mutation_rand_mats2 = h5read("/data/dsge_data_dir/dsgejl/interns2017/mutation_rand_mats.h5",
                                 "mutation_rand_mats2")
    s_lag_tempered_rand_mat = h5read("/data/dsge_data_dir/dsgejl/interns2017/s_lag_tempered.h5",
                                "s_lag_tempered_rand_mat")
    ε_rand_mats = h5read("/data/dsge_data_dir/dsgejl/interns2017/eps_rand_mats.h5",
                                "eps_rand_mats")
    
    if deterministic
        # Testing: Generate pre-defined random shocks for s and ε
        # s_lag_tempered_rand_mat = randn(n_states, n_particles)
        #ε_rand_mat = ε_rand_mats[:,:,eps_id] #randn(n_observables, n_particles)
        #path = dirname(@__FILE__)  
        
        #h5open("$path/../../../test/reference/matricesForTPF.h5", "w") do file
        #    write(file, "s_lag_tempered_rand_mat", s_lag_tempered_rand_mat)
        #    write(file, "eps_rand_mat", ε_rand_mat)
        #end
    else
        # Draw initial particles from the distribution of s₀: N(s₀, P₀) 
        s_lag_tempered_rand_mat = randn(n_states, n_particles)
    end

    ### Change back to get_chol later!!
    s_lag_tempered = repmat(s0, 1, n_particles) + Matrix(chol(P0))'*s_lag_tempered_rand_mat
    @show chol(P0)
    @show s0
    @show sum(s_lag_tempered,2)


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
        nonmissing          = !isnan(y_t)
        y_t                 = y_t[nonmissing]        
        ZZ_t                = ZZ[nonmissing,:]
        DD_t                = DD[nonmissing]
        EE_t                = EE[nonmissing,nonmissing]
        QQ_t                = QQ[nonmissing,nonmissing]
        RRR_t               = RRR[:,nonmissing]
        sqrtS2_t            = RRR_t*get_chol(QQ_t)'
        n_observables_t     = length(y_t)
        ε                   = zeros(n_observables_t)
        mutation_rand_mat_t = mutation_rand_mat[nonmissing,:] 

        if deterministic # When not testing, want a new random ε each time 
            ε_rand_mat = ε_rand_mats[:,:,T]
        else    
            ε_rand_mat = randn(n_observables_t, n_particles)
        end

        # TEMPORARY
        #ε_rand_mat = randn(n_observables_t, n_particles)

        # Forecast forward one time step
        @show sum(sum(s_lag_tempered))
        s_t_nontempered = TTT*s_lag_tempered + sqrtS2_t*ε_rand_mat
        @show sum(sum(s_t_nontempered))

        # Error for each particle
        p_error = repmat(y_t - DD_t, 1, n_particles) - ZZ_t*s_t_nontempered
        @show sum(sum(p_error))

        # Solve for initial tempering parameter φ_1
        if deterministic 
            φ_1 = 0.25
        else
            init_Ineff_func(φ) = solve_inefficiency(φ, 2.0*pi, y_t, p_error, EE_t, 
                                                    initialize=true)-rstar
            φ_1 = fzero(init_Ineff_func, 1e-7, 1.0, xtol=xtol)
        end
        
         if VERBOSITY[verbose] >= VERBOSITY[:low]
            @show φ_1 
            println("------------------------------")
        end
              
        # Update weights array and resample particles
        loglik, weights, s_lag_tempered, ε, id = correct_and_resample2!(φ_1, 0.0, y_t, p_error,
                                                     incremental_weights, weights, s_lag_tempered,
                                                     ε_rand_mat, EE_t, n_particles, deterministic,
                                                     initialize=true,id_cur=id_cur,matlab_resampled_ids=matlab_resampled_ids)

#        @show loglik
#        @show weights
#        @show s_lag_tempered
#        @show ε
#        @show id

        #resampling_ids[ids_i,:] = id
        #ids_i += 1
        id_cur += 1 
        # Update likelihood
        lik[t] += loglik
        @show lik[t]
        # Tempering Initialization
        count = 2 # Accounts for initialization and final mutation
        φ_old = φ_1

        # First propagation
        s_t_nontempered = TTT*s_lag_tempered + sqrtS2_t*ε
        p_error = repmat(y_t - DD_t, 1, n_particles) - ZZ_t*s_t_nontempered         
        
        if deterministic
            ineff_check = rstar + 1
        else
            ineff_check = solve_inefficiency(1.0, φ_1, y_t, p_error, EE_t)         
        end

        if VERBOSITY[verbose] >= VERBOSITY[:high]
            deterministic ? println("Deterministic:") : println("Non-deterministic:")
            @show ineff_check
        end

        #--------------------------------------------------------------
        # Main Algorithm
        #--------------------------------------------------------------
        while ineff_check > rstar

            # Define inefficiency function
            init_ineff_func(φ) = solve_inefficiency(φ, φ_old, y_t, p_error, EE_t) - rstar
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
                    ineff_check = solve_inefficiency(1.0, φ_new, y_t, p_error, EE_t)
                end
               
                if VERBOSITY[verbose] >= VERBOSITY[:low]
                    @show φ_new
                end

                # Update weights array and resample particles
                loglik, weights, s_lag_tempered, ε, id = correct_and_resample2!(φ_new, φ_old, y_t,
                                                             p_error, incremental_weights, weights,
                                                             s_lag_tempered, ε, EE_t, n_particles,
                                                             deterministic, id_cur=id_cur,matlab_resampled_ids=matlab_resampled_ids)
                #resampling_ids[ids_i,:] = id
                #ids_i += 1
                id_cur += 1
                # Update likelihood
                lik[t] += loglik
                @show lik[t]
                # Update value for c
                c = update_c!(m, c, accept_rate, target)
                
                if VERBOSITY[verbose] >= VERBOSITY[:low]
                    if VERBOSITY[verbose] >= VERBOSITY[:high]
                        @show c
                    end
                    println("------------------------------")
                end
                                
                # Mutation Step
                accept_vec = zeros(n_particles)
                print("Mutation ")        
                tic()

                if parallel
                    print("(in parallel) ")                    
                    #out = pmap(i -> mutation(system, y_t, s_lag_tempered[:,i], ε[:,i], c, 
                    #           N_MH, deterministic, nonmissing, mutation_rand_mat), 1:n_particles)
                    out = @sync @parallel (hcat) for i=1:n_particles
                    mutation_matlab(system, y_t, s_lag_tempered[:,i], ε[:,i], c, N_MH, deterministic, nonmissing, mutation_rand_mat_t, mutation_rand_mats1[:,id_m,i], mutation_rand_mats2[:,id_m,i])
                    end
                else 
                    print("(not parallel) ")
                    out = [mutation(system, y_t, s_lag_tempered[:,i], ε[:,i], c, N_MH, 
                                    deterministic, nonmissing, mutation_rand_mat_t, mutation_rand_mats1[:,id_m,i], mutation_rand_mats2[:,id_m,i]) 
                           for i=1:n_particles]
                end
                times[t] = toc()                
                id_m += 1

                for i = 1:n_particles
                    s_t_nontempered[:,i] = out[i][1]
                    ε[:,i] = out[i][2]
                    accept_vec[i] = out[i][3]
                end

                # Calculate average acceptance rate
                accept_rate = mean(accept_vec)

                # Get error for all particles
                p_error = repmat(y_t-DD_t, 1, n_particles)-ZZ_t*s_t_nontempered
                φ_old = φ_new
                len_phis[t] += 1

            # If no solution exists within interval, set inefficiency to rstar
            else 
                if VERBOSITY[verbose] >= VERBOSITY[:high]
                    println("No solution in interval.")
                end
                ineff_check = rstar
            end
            gc()

            # With phi schedule, exit after one iteration, thus set ineff_check=0
            if deterministic
                ineff_check = 0.0
            end

            if VERBOSITY[verbose] >= VERBOSITY[:high]
                @show ineff_check
            end
        end

        if VERBOSITY[verbose] >= VERBOSITY[:high]
            println("Out of main while-loop.")
        end
        
        #--------------------------------------------------------------
        # Last Stage of Algorithm: φ_new := 1.0
        #--------------------------------------------------------------
        φ_new = 1.0

        # Update weights array and resample particles.
        loglik, weights, s_lag_tempered, ε, id = correct_and_resample2!(φ_new, φ_old, y_t, p_error,
                                                     incremental_weights, weights, s_lag_tempered,
                                                     ε, EE_t, n_particles, deterministic, id_cur=id_cur,matlab_resampled_ids=matlab_resampled_ids)
        id_cur += 1
        #resampling_ids[ids_i,:] = id
        #ids_i += 1

        # Update likelihood
        lik[t] += loglik
        @show lik[t]
        # Update c
        c = update_c!(m, c, accept_rate, target)
        
        # Final round of mutation
        accept_vec = zeros(n_particles)

        if parallel
            # out = pmap(i -> mutation(system, y_t, s_lag_tempered[:,i], ε[:,i], c, N_MH, 
            #            deterministic,nonmissing,mutation_rand_mat), 1:n_particles)
            out = @sync @parallel (hcat) for i=1:n_particles
                mutation_matlab(system, y_t, s_lag_tempered[:,i], ε[:,i], c,N_MH, deterministic,
                         nonmissing, mutation_rand_mat_t, mutation_rand_mats1[:,id_m,i], mutation_rand_mats2[:,id_m,i])
            end
        else 
            out = [mutation(system, y_t, s_lag_tempered[:,i], ε[:,i], c, N_MH, deterministic,
                            nonmissing, mutation_rand_mat_t, mutation_rand_mats1[:,id_m,i], mutation_rand_mats2[:,id_m,i]) for i=1:n_particles]
        end
        
        id_m += 1
        # Unwrap output
        for i = 1:n_particles
            s_t_nontempered[:,i] = out[i][1]
            ε[:,i] = out[i][2]
            accept_vec[i] = out[i][3]
        end
        
        # Store for next time iteration
        accept_rate = mean(accept_vec)

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
update_c!(m::AbstractModel, c_in::Float64, accept_in::Float64, target_in::Float64)
```
Update value of c by expression that is function of the target and mean acceptance rates.
Returns the new c, in addition to storing it in the model settings.

"""
function update_c!(m::AbstractModel,c_in::Float64, accept_in::Float64, target_in::Float64)
    c_out = c_in*(0.95 + 0.1*exp(20*(accept_in - target_in))/(1 + exp(20*(accept_in - target_in))))
    return c_out
end

"""
```
correct_and_resample!(φ_new::Float64, φ_old::Float64, y_t::Array, p_error::Array,incremental_weights::Array,weights::Array, s_lag_tempered::Array, ε::Array, EE::Array, n_particles::Int64; initialize::Bool=false)
```
Calculate densities, normalize and reset weights, call multinomial resampling, update state and error vectors,reset error vectors to 1,and calculate new log likelihood.
Returns log likelihood, weight, state, and ε vectors.

"""
function correct_and_resample2!(φ_new::Float64, φ_old::Float64, y_t::Array{Float64,1}, p_error::Array{Float64,2}, incremental_weights::Array{Float64,1}, weights::Array{Float64,1}, s_lag_tempered::Array{Float64,2}, ε::Array{Float64,2}, EE::Array{Float64,2}, n_particles::Int64, deterministic::Bool; initialize::Bool=false,id_cur::Int=1,matlab_resampled_ids=zeros(1,1))
    # Calculate initial weights
    for n=1:n_particles
        incremental_weights[n] = incremental_weight(φ_new, φ_old, y_t, p_error[:,n], EE, 
                                                    initialize=initialize)
    end   

    # Normalize weights
    weights = (incremental_weights.*weights) ./ mean(incremental_weights.*weights)
    
    # Resampling
    if deterministic
        # id = seeded_multinomial_resampling(weights)
        id = matlab_resampled_ids[:,id_cur]
    else
        id = multinomial_resampling(weights)
    end

    # Update arrays for resampled indices
    s_lag_tempered = s_lag_tempered[:,id]
    ε = ε[:,id]

    # Reset weights to ones
    weights = ones(n_particles)

    # Calculate likelihood
    loglik = log(mean(incremental_weights.*weights))
    
    return loglik, weights, s_lag_tempered, ε, id
end

"""
```
incremental_weight{S<:Float64}(φ_new::S, φ_old::S, y_t::Array{S,1}, p_error::Array{S,1}, 
    EE::Array{S,2}; initialize::Bool=false)
```
### Arguments
-`φ_new::S`: 
-`φ_old::S`: φ value before last
-`y_t::Array{S,1}`: Vector of observables for time t
-`p_error::Array{S,1}`: A single particle's error: y_t - Ψ(s_t)
-`EE::Array{S,2}`: Measurement error covariance matrix
-`initialize::Bool`: Flag indicating whether one is solving for incremental weights during 
                     the initialization of weights; default is false.

### Output

Returns the probability evaluated at p_error.
"""
function incremental_weight{S<:Float64}(φ_new::S, φ_old::S, y_t::Array{S,1}, p_error::Array{S,1}, 
                                   EE::Array{S,2}; initialize::Bool=false)

    # Initialization step (using 2π instead of φ_old)
    if initialize
        return (φ_new/(2*pi))^(length(y_t)/2) * (det(EE)^(-1/2)) * 
            exp(-1/2 * p_error' * φ_new * inv(EE) * p_error)[1]
    
    # Non-initialization step (tempering and final iteration)
    else
        return (φ_new/φ_old)^(length(y_t)/2) * 
            exp(-1/2 * p_error' * (φ_new - φ_old) * inv(EE) * p_error)[1]
    end
end
