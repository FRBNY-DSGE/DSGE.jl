using DSGE
function tpf_loop{S<:AbstractFloat}(m::AbstractModel, yy::Array, system::System{S},
    s0::Array{S}, P0::Array; verbose::Symbol=:low, include_presample::Bool=true)

    # Unpack system
    RRR = system[:RRR]
    TTT = system[:TTT]
    EE  = system[:EE]
    DD  = system[:DD]
    ZZ  = system[:ZZ]
    QQ  = system[:QQ]    
    
    # Get tuning parameters from the model
    r_star       = get_setting(m, :tpf_r_star)
    c            = get_setting(m, :tpf_c_star)
    accept_rate  = get_setting(m, :tpf_accept_rate)
    target       = get_setting(m, :tpf_target)
    N_MH         = get_setting(m, :tpf_n_mh_simulations)
    n_particles  = get_setting(m, :tpf_n_particles)
    adaptive     = get_setting(m, :tpf_adaptive)
    xtol         = get_setting(m, :tpf_x_tolerance)
    parallel     = get_setting(m, :use_parallel_workers)

    # Determine presampling periods
    n_presample_periods = (include_presample) ? 0 : get_setting(m, :n_presample_periods)
    
    # End time (last period)
    T = size(yy,2) 
    
    # Initialization
    n_observables       = size(QQ,1)
    n_states            = size(ZZ,2)
    lik                 = zeros(T)
    Neff                = zeros(T)
    n_φ_steps           = ones(T)
    times               = zeros(T)
    weights             = ones(n_particles)
    incremental_weights = zeros(n_particles)
    times               = zeros(T)
        
    # Draw initial particles from the distribution of s₀: N(s₀, P₀) 
    s_lag_tempered = repmat(s0, 1, n_particles) + Matrix(chol(P0))'*randn(n_states, n_particles)
   
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
        nonmissing      = !isnan(y_t)
        y_t             = y_t[nonmissing]        
        ZZ_t            = ZZ[nonmissing,:]
        DD_t            = DD[nonmissing]
        EE_t            = EE[nonmissing,nonmissing]
        QQ_t            = QQ[nonmissing,nonmissing]
        RRR_t           = RRR[:,nonmissing]
        sqrtS2_t        = RRR_t*get_chol(QQ_t)'
        n_observables_t = length(y_t)
        #ε               = zeros(n_observables_t)
        
        # Draw random shock ε
        ε = randn(n_observables_t, n_particles)

        # Forecast forward one time step
        s_t_nontempered = TTT*s_lag_tempered + sqrtS2_t*ε

        converged = false
        φ_old = 0.001
        φ_min = 0.001
        
        oldweight = ones(n_particles)
        while !converged
            tic()
            # Error for each particle
            p_error = repmat(y_t - DD_t, 1, n_particles) - ZZ_t*s_t_nontempered
            #@show mean(p_error,2)
   #=         incwt=ones(n_particles)
            
            for i=1:n_particles
                incwt[i] = -.5*(p_error[:,i]'*inv(EE)*p_error[:,i])[1]
                #check what's going on here
                oldweight[i] = weights[i]
            end
            println("old weights")
            @show mean(oldweight)
            wt_kernel = incwt
            ess0 = ineff_func_h(φ_min, φ_old, wt_kernel, n_particles,r_star)
            ess1 = ineff_func_h(1.0, φ_old, wt_kernel, n_particles,r_star)
            
            @show mean(wt_kernel)
            if ess1 < 0.0
                φ_new = 1.0
                converged=true
            elseif ess0>0
                print("temp high")
            else
                func(φ) = ineff_func_h(φ, φ_old,wt_kernel,n_particles,r_star)
                φ_new = fzero(func, φ_min, 1.0, xtol=0.1)
            end
            for i = 1:n_particles
                if φ_old>0.0 
                    weights[i] = oldweight[i]*(φ_new/φ_old)^(length(y_t)/2)*exp(wt_kernel[i]*(φ_new-φ_old))
                else
                    weights[i] = (2*pi)^(length(y_t)/2)*(det(EE))^(-1/2)*φ_new^(length(y_t)/2)*exp(φ_new*wt_kernel[i])*oldweight[i]
                end
            end
=#
            initialize = (φ_old==0.001)
     #=       @show initialize
            for j=1:n_particles
                weights[j] = incremental_weight(φ_new, φ_old,y_t,p_error[:,j],EE_t,initialize=initialize)
            end
            =#          
            init_Ineff_func(φ) = solve_inefficiency(φ,φ_old,y_t,p_error,EE_t,initialize=initialize)-r_star
            if prod(sign([init_Ineff_func(φ_old),init_Ineff_func(1.0)])) ==-1
                φ_new = fzero(init_Ineff_func,φ_old,1.0,xtol=xtol)
                ineff_check = solve_inefficiency(1.0,φ_old,y_t,p_error,EE_t)
                if ineff_check <=r_star
                    φ_new = 1.0
                    converged=true
                end
                for j=1:n_particles
                    weights[j]=incremental_weight(φ_new,φ_old,y_t,p_error[:,j],EE,initialize=initialize)
                end

                println("before normalize")
                @show mean(weights)
                weights = weights./mean(weights)        
                id = multinomial_resampling(weights)
                s_lag_tempered = s_lag_tempered[:,id]
                ε= ε[:,id]
               @show mean(weights)
                loglik= log(mean(weights))
                weights=ones(n_particles)
         
                φ_old = φ_new
                φ_min = φ_old
        
                lik[t] += loglik
               
                c= update_c!(c,accept_rate,target)
                accept_vec = zeros(n_particles)
                out = @sync @parallel (hcat) for i=1:n_particles
                    mutation(system,y_t,s_lag_tempered[:,i],ε[:,i],c,N_MH,nonmissing) 
                end
                for i=1:n_particles
                    s_t_nontempered[:,i] = out[i][1]
                    ε[:,i]=out[i][2]
                    accept_vec[i]=out[i][3]
                end
                accept_rate = mean(accept_vec)
                φ_old = φ_new
                
                #mutate
        #=        s_init = s_lag_tempered
                @show φ_new
                for k = 1: N_MH
                    ε_new = ε + c*randn(length(y_t),n_particles)
                    acpt = 0
                    for j = 1:n_particles
                        s_new_e = TTT*s_init[:,j]+sqrtS2_t*ε_new[:,j]
                        s_init_e = TTT*s_init[:,j] + sqrtS2_t*ε[:,j]
                        error_new = y_t - ZZ_t*s_new_e - DD_t
                        error_init = y_t - ZZ_t*s_init_e - DD_t
                        
                        post_new = log(pdf(MvNormal(zeros(length(y_t)),EE_t),error_new)*pdf(MvNormal(zeros(length(ε_new[:,j])),eye(length(ε_new[:,j]))),ε_new[:,j]))
                        post_init = log(pdf(MvNormal(zeros(length(y_t)),EE_t),error_init)[1]*pdf(MvNormal(zeros(length(ε[:,j])),eye(length(ε[:,j]))),ε[:,j])[1])
                        α = exp(φ_new*(post_new-post_init))
                        if rand() <α
                            ε = ε_new
                            s_t_nontempered[:,j] = s_new_e
                            acpt += 1
                        end
                    end
                    acpt_rate = sum(acpt)/n_particles
                    c = c*(0.95+0.1*exp(20*(acpt_rate-0.4))/(1+exp(20*(acpt_rate-0.4))))
                end
=#               
            else 
                converged=true
            end
        end
            # Update value for c
            # c = update_c!(c, accept_rate, target)
            
            # Mutation Step
            #accept_vec = zeros(n_particles)
            #print("Mutation ")        
            # tic()
            
            # print("(in parallel) ")                    
            # out = @sync @parallel (hcat) for i=1:n_particles
            #      mutation_loop(system, y_t, s_lag_tempered[:,i], ε[:,i], c, N_MH, nonmissing)
            ## end
            #toc()                
            #for i = 1:n_particles
            #   s_t_nontempered[:,i] = out[i][1]
            #  ε[:,i] = out[i][2]
            #  accept_vec[i] = out[i][3]
            # end            
    
        n_φ_steps[t] += 1
            
    Neff[t] = (n_particles^2)/sum(weights.^2)

    s_lag_tempered = s_t_nontempered
    print("Completion of one period ") 
    toc()
    @show lik[t]

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("=============================================")
    end
end
    # Return vector of likelihood indexed by time step and Neff
    return Neff[n_presample_periods + 1:end], lik[n_presample_periods + 1:end]
end




"""
```
get_chol(mat::Aray)
```
Calculate and return the Cholesky of a matrix.

"""
@inline function get_chol(mat::Array{Float64,2})
    return Matrix(chol(nearestSPD(mat)))
end

"""
```
update_c!(m::AbstractModel, c_in::Float64, accept_in::Float64, target_in::Float64)
```
Update value of c by expression that is function of the target and mean acceptance rates.
Returns the new c, in addition to storing it in the model settings.

"""
@inline function update_c!{S<:Float64}(c_in::S, accept_in::S, target_in::S)
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
function correct_and_resample!{S<:Float64}(φ_new::S, φ_old::S, y_t::Array{S,1}, p_error::Array{S,2}, incremental_weights::Array{S,1}, weights::Array{S,1}, s_lag_tempered::Array{S,2}, ε::Array{S,2}, EE::Array{S,2}, n_particles::Int; initialize::Bool=false)
    
    # Calculate initial weights
    for n=1:n_particles
        incremental_weights[n] = incremental_weight(φ_new, φ_old, y_t, p_error[:,n], EE, 
                                                    initialize=initialize)
    end   

    # Normalize weights
    weights = (incremental_weights.*weights) ./ mean(incremental_weights.*weights)
    
    # Resampling
    id = multinomial_resampling(weights)
    
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

function ineff_func_h(φ,φ_old,incwt,n_particles,rstar)
    if φ_old==0 
        new_weight = exp(φ*incwt)
        new_weight2 = exp(2*φ*incwt)
    else 
        new_weight = exp((φ-φ_old)*incwt)
        new_weight2= exp(2*(φ-φ_old)*incwt)
    end
    a1 = sum(new_weight2)/n_particles
    a2 = (sum(new_weight)/n_particles)^2
    return (a1/a2)-rstar
end        