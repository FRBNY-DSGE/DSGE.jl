using ClusterManagers, HDF5
using Roots: fzero
using QuantEcon: solve_discrete_lyapunov

addprocs_sge(30, queue="background.q")

@everywhere using DSGE

function the_problem()

m = SmetsWouters("ss1", testing=true)

path = dirname(@__FILE__)
filesw= "/data/dsge_data_dir/dsgejl/realtime/input_data/data"
data = readcsv("$filesw/realtime_spec=smets_wouters_hp=true_vint=110110.csv",header=true)
data = convert(Array{Float64,2},data[1][:,2:end])
data=data'
params = h5read("$filesw/../../output_data/smets_wouters/ss0/estimate/raw/paramsmode_vint=110110.h5","params")
 push!(params, m[:e_y].value, m[:e_L].value, m[:e_w].value, m[:e_π].value, m[:e_R].value, m[:e_c].value, m[:e_i].value)
update!(m,params)
system = compute_system(m)

mats = "$path/../../test/reference/sysMats.h5"
h5open(mats, "w") do file
    write(file, "R", system.transition.RRR)
    write(file, "T", system.transition.TTT)
    write(file, "E", system.measurement.EE)
    write(file, "D", system.measurement.DD)
    write(file, "Z", system.measurement.ZZ)
    write(file, "Q", system.measurement.QQ)
    write(file, "sqrtS2",system.transition.RRR*get_chol(system.measurement.QQ)')
end

RRR = h5read(mats, "R")
TTT = h5read(mats, "T")
EE = h5read(mats, "E")
DD = h5read(mats, "D")
ZZ = h5read(mats, "Z")
S2 = h5read(mats, "Q")
QQ = h5read(mats, "Q")
sqrtS2 = h5read(mats, "sqrtS2")


rand_mat = randn(size(S2,1),1)
m<=Setting(:tpf_rand_mat,rand_mat)

m<=Setting(:tpf_rstar,2.0)
m<=Setting(:tpf_N_MH, 10)
m<=Setting(:tpf_c,0.1)
m<=Setting(:tpf_acpt_rate,0.5)
m<=Setting(:tpf_trgt,0.25)
n_particles=4000
m<=Setting(:tpf_n_particles,n_particles)
m<=Setting(:use_parallel_workers, true)
m<=Setting(:x_tolerance, zero(float(0)))
m<=Setting(:tpf_deterministic, false)

s0 = zeros(size(TTT)[1])
P0= nearestSPD(solve_discrete_lyapunov(TTT,RRR*S2*RRR'))

#tpf(m,data,system,s0,P0)

# Store model parameters
    
    # Get tuning parameters from the model
    rstar         = get_setting(m,:tpf_rstar)
    c             = get_setting(m,:tpf_c)
    acpt_rate     = get_setting(m,:tpf_acpt_rate)
    trgt          = get_setting(m,:tpf_trgt)
    N_MH          = get_setting(m,:tpf_n_mh_simulations)
    n_particles   = get_setting(m,:tpf_n_particles)
    deterministic = get_setting(m,:tpf_deterministic)
    xtol          = get_setting(m,:tpf_x_tolerance) # Tolerance of fzero
    parallel      = get_setting(m,:use_parallel_workers) # Get setting of parallelization

# End time (last period)
    T = size(data,2)
    # Size of covariance matrix
    n_errors = size(QQ,1)
    # Number of states
    n_states = size(ZZ,2)
    
    # Initialization
    lik         = zeros(T)
    Neff        = zeros(T)
    len_phis    = ones(T)
    weights     = ones(n_particles)
    density_arr = zeros(n_particles)

s_lag_tempered_rand_mat = randn(n_states,n_particles)
        ε_rand_mat = randn(n_errors, n_particles)

s_lag_tempered = repmat(s0,1,n_particles) + get_chol(P0)'*s_lag_tempered_rand_mat

for t=1:T
    @show t
    tic()
    yt = data[:,25]
    nonmissing =!isnan(yt)

    s_t_nontempered = TTT*s_lag_tempered + sqrtS2*ε_rand_mat
        
    # Error for each particle
    perror = repmat(yt-DD,1,n_particles)-ZZ*s_t_nontempered

    init_Ineff_func(φ) = ineff_func(φ, 2.0*pi, yt, perror, EE, initialize=1)-rstar
    φ_1 = fzero(init_Ineff_func,0.00000001, 1.0, xtol=xtol)

    loglik, weights, s_lag_tempered, ε = correct_and_resample(φ_1,0.0,yt,perror,density_arr,weights,s_lag_tempered,ε_rand_mat,EE,n_particles,deterministic,initialize=1)
        
    # Update likelihood
    lik[t]=lik[t]+loglik

        # Tempering Initialization
        count = 2 # Accounts for initialization and final mutation
        φ_old = φ_1

        # First propagation
        s_t_nontempered = TTT*s_lag_tempered + sqrtS2*ε
        perror = repmat(yt-DD, 1, n_particles) - ZZ*s_t_nontempered         
    check_ineff=ineff_func(1.0, φ_1, yt, perror, EE)         

    while (check_ineff>rstar)

        # Define inefficiency function
        init_ineff_func(φ) = ineff_func(φ, φ_old, yt, perror, EE)-rstar
        φ_interval = [φ_old, 1.0]
        fphi_interval = [init_ineff_func(φ_old) init_ineff_func(1.0)]

        count += 1

        # Check solution exists within interval
        if prod(sign(fphi_interval))== -1 || deterministic
                
           
            # Set φ_new to the solution of the inefficiency function over interval
            φ_new = fzero(init_ineff_func, φ_interval, xtol=xtol)
            check_ineff = ineff_func(1.0, φ_new, yt, perror, EE, initialize=0)
            
            # Update weights array and resample particles
            loglik, weights, s_lag_tempered, ε = correct_and_resample(φ_new,φ_old,yt,perror,density_arr,weights,s_lag_tempered,ε,EE,n_particles, deterministic, initialize=0)
                
            # Update likelihood
            lik[t] = lik[t] + loglik
                
            # Update value for c
            c = update_c(m,c,acpt_rate,trgt)
                                            
            # Update covariance matrix
            μ = mean(ε,2)
            cov_s = (1/n_particles)*(ε-repmat(μ,1,n_particles))*(ε - repmat(μ,1,n_particles))'
              
            if !isposdef(cov_s)
                cov_s = diagm(diag(cov_s))
            end
            
            # Mutation Step
            acpt_vec=zeros(n_particles)
            print("Mutation ")
                
            tic()

            if parallel
                print("(in parallel) ")
                c = get_setting(m,:tpf_c)
                N_MH=get_setting(m,:tpf_n_mh_simulations)
                deterministic=get_setting(m,:tpf_deterministic)
                #out = pmap(i->mutation(c, N_MH,deterministic,yt,s_lag_tempered[:,i],ε[:,i],cov_s,nonmissing), 1:n_particles)
                out = @sync @parallel (hcat) for i=1:n_particles
                    mutation(c,N_MH,deterministic,system,yt,s_lag_tempered[:,i],ε[:,i],cov_s,nonmissing)
                end
            else
                print("(not parallel)")
                out = [mutation(c,N_MH,deterministic,system,yt,s_lag_tempered[:,i],ε[:,i],cov_s,nonmissing) for i=1:n_particles]
            end               
            toc()
            for i = 1:n_particles
                s_t_nontempered[:,i] = out[i][1]
                ε[:,i] = out[i][2]
                acpt_vec[i]=out[i][3]
            end
            acpt_rate = mean(acpt_vec)

            # Get error for all particles
            perror = repmat(yt-DD, 1, n_particles)-ZZ*s_t_nontempered
            φ_old = φ_new
            len_phis[t]+=1

            # If no solution exists within interval, set inefficiency to rstar
        else 
            check_ineff = rstar
        end
    end

            
    #————————————————————–
    # Last Stage of Algorithm: φ_new=1
    #————————————————————–
    φ_new = 1.0

    # Update weights array and resample particles.
    loglik, weights, s_lag_tempered, ε = correct_and_resample(φ_new,φ_old,yt,perror,density_arr,weights,s_lag_tempered,ε,EE,n_particles,deterministic,initialize=0)

    # Update likelihood
    lik[t]=lik[t]+loglik

    c = update_c(m,c,acpt_rate,trgt)
    μ = mean(ε,2)
    
    cov_s = (1/n_particles)*(ε-repmat(μ,1,n_particles))*(ε - repmat(μ,1,n_particles))'
                
    if isposdef(cov_s)== false 
        cov_s = diagm(diag(cov_s))
    end
        
    # Final round of mutation
    acpt_vec = zeros(n_particles)
    
    if parallel
        c = get_setting(m,:tpf_c)
        N_MH=get_setting(m,:tpf_n_mh_simulations)
        deterministic=get_setting(m,:tpf_deterministic)
        #out = pmap(i -> mutation(c,N_MH,deterministic,system,yt,s_lag_tempered[:,i],ε[:,i],cov_s,nonmissing), 1:n_particles)
        out = @sync @parallel (hcat) for i=1:n_particles
            mutation(c,N_MH,deterministic,system,yt,s_lag_tempered[:,i],ε[:,i],cov_s,nonmissing)
        end
    else 
        out = [mutation(c,N_MH,deterministic,system,yt,s_lag_tempered[:,i],ε[:,i],cov_s,nonmissing) for i=1:n_particles]
    end

    for i = 1:n_particles
        s_t_nontempered[:,i] = out[i][1]
        ε[:,i] = out[i][2]
        acpt_vec[i] = out[i][3]
    end
        
    # Store for next time iteration
    acpt_rate = mean(acpt_vec)

    Neff[t] = (n_particles^2)/sum(weights.^2 )
    s_lag_tempered = s_t_nontempered
    print("Completion of one period ")
    gc()
    toc()
end
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
    c_out = c_in*(0.95+0.1*exp(16*(acpt_in-trgt_in))/(1+exp(16*(acpt_in-trgt_in))))
    m <= Setting(:tpf_c, c_out)
    return c_out
end

"""
```
correct_and_resample(φ_new::Float64, φ_old::Float64, yt::Array, perror::Array,density_arr::Array,weights::Array, s_lag_tempered::Array, ε::Array, EE::Array, n_particles::Int64; initialize::Int64=0)
```
Calculate densities, normalize and reset weights, call multinomial resampling, update state and error vectors,reset error vectors to 1,and calculate new log likelihood.
Returns log likelihood, weight, state, and ε vectors.
"""
function correct_and_resample(φ_new::Float64, φ_old::Float64, yt::Array{Float64,1}, perror::Array{Float64,2}, density_arr::Array{Float64,1}, weights::Array{Float64,1}, s_lag_tempered::Array{Float64,2}, ε::Array{Float64,2}, EE::Array{Float64,2}, n_particles::Int64, deterministic::Bool; initialize::Int64=0)

    # Calculate initial weights
    for n=1:n_particles
        density_arr[n]=density(φ_new, φ_old, yt, perror[:,n], EE, initialize=initialize)
    end   

    # Normalize weights
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
    
    return loglik, weights, s_lag_tempered, ε
end

function zlb_regime_indices{S<:AbstractFloat}(m::AbstractModel{S},data::Matrix{S})
    #make sure the data matrix has all time periods when passing in or this won't work
    T = size(data,2)
    if n_anticipated_shocks(m) >0
        regime_inds = Vec{Range{Int64}}(2)
        regime_inds[1] = 1:index_zlb_start(m)-1
        regime_inds[2] = index_zlb_start(m):T 
    else 
        regime_inds = Range{Int64}[1:T]
    end
end

function zlb_regime_matrices{S<:AbstractFloat}(m::AbstractModel{S},system::System{S})
    if !all(x -> x==0, system[:MM])
        error("previously this error said Kalman filter and smoothers not implemented for nonzero MM however i'm not sure if this still applies to the TPF")
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
    ZZs = fill(system[:ZZ], n_regimes)
    DDs = fill(system[:DD], n_regimes)
    EEs = fill(system[:EE], n_regimes)

    return TTTs, RRRs, CCCs, QQs, ZZs, DDs, EEs
end