"""
```
smc(m::AbstractModel,data::Matrix)
smc(m::AbstractModel,data::DataFrame)
smc(m::AbstractModel)
```

### Arguments:

- `m`: A model object, from which its parameters values, prior dists, and bounds will be referenced
- `data`: A matrix or data frame containing time series of the observables to be used in the calculation of the posterior/likelihood

### Keyword Arguments:
- `verbose`: The desired frequency of function progress messages printed to standard out.
	- `:none`: No status updates will be reported.
	- `:low`: Status updates for SMC initialization and recursion will be included.
	    - `:high`: Status updates for every iteration of SMC is output, which includes the mu and sigma of each individual parameter draw after each iteration, as well as the calculated acceptance rate, the ESS, and the number of times resampled.

### Outputs

- `mu`: The mean of a particular parameter as calculated by taking the weighted sum of that parameter from the particle draws made in that particular stage Nφ.
- `sig`: The standard deviation of a particular parameter.

### Overview

Sequential Monte Carlo can be used in lieu of Random Walk Metropolis Hastings to generate parameter samples from high-dimensional parameter spaces using sequentially constructed proposal densities.

The implementation here is based upon Edward Herbst and Frank Schorfheide's 2014 paper 'Sequential Monte Carlo Sampling for DSGE Models' and the code accompanying their book 'Bayesian Estimation of DSGE Models'.

SMC is broken up into three main steps:

- `Correction`: Reweight the particles from stage n-1 by defining "incremental weights", incweight, which gradually incorporate the likelihood function p(Y|θ(i,n-1)) into the particle weights.
- `Selection`: Resample the particles if the distribution of particles begins to degenerate, according to a tolerance level ESS < N/2. The current resampling technique employed is systematic resampling.

- `Mutation`: Propagate particles {θ(i),W(n)} via N(MH) steps of a Metropolis Hastings algorithm.
"""
function smc(m::AbstractModel, data::Matrix; verbose::Symbol=:low)
    #--------------------------------------------------------------
    #Set Parameters of Algorithm
    #--------------------------------------------------------------

    n_Φ = get_setting(m, :n_Φ)
    n_part = get_setting(m, :n_particles)
    n_params = n_parameters(m)
    n_blocks = get_setting(m, :n_smc_blocks)

    fixed_para_inds = find([ θ.fixed for θ in m.parameters])
    free_para_inds  = find([!θ.fixed for θ in m.parameters])
    c = get_setting(m, :step_size_smc)
    accept = get_setting(m, :init_accept)
    target = get_setting(m, :target_accept)
    
    draw_type = get_setting(m,:draw_type)
    resampler(x) = if get_setting(m, :resampler_smc) == :systematic
        systematic_resampling(m, x)
    elseif get_setting(m, :resampler_smc) == :multinomial
        sample(1:n_part, weights(x), n_part, replace = true)
    else
        throw("Invalid resampler in SMC. Set :resampler_smc setting to either :systematic or :multinomial")
    end

    parallel = get_setting(m, :use_parallel_workers)

    # Creating the tempering schedule
    λ = get_setting(m, :λ)
    tempering_schedule = ((collect(1:1:n_Φ)-1)/(n_Φ-1)).^λ

    #Matrices for storing

    para_sim = zeros(n_Φ, n_part, n_params) # parameter draws
    weight_sim = zeros(n_part, n_Φ) # weights
    zhat = zeros(n_Φ, 1) # normalization constant
    n_resamp = 0 # record # of iteration resampled

    c_sim = zeros(n_Φ,1) # scale parameter
    ESS_sim = zeros(n_Φ,1) # ESS
    accept_sim = zeros(n_Φ,1) # average acceptance rate
    rsmp_sim = zeros(n_Φ,1) # 1 if resampled

    #--------------------------------------------------------------
    #Initialize Algorithm: Draws from prior
    #--------------------------------------------------------------

    if VERBOSITY[verbose] >= VERBOSITY[:low]
	println("\n\n SMC starts ....  \n\n  ")
    end

    # Particle draws from the parameter's marginal priors
    prior_sim = zeros(n_part,n_params)

    # Posterior values at prior draws
    loglh = zeros(n_part)
    logpost = zeros(n_part)

    # fix seed if testing
    if m.testing
        println("seed set")
        srand(42)
    end

    ##Get parameters at mode and hessian
    params = h5open(rawpath(m,"estimate", "paramsmode.h5"),"r") do file
        read(file, "params")
    end

    # fn = hessian_path(m)
    fn = rawpath(m,"estimate","hessian.h5")
    hessian = h5open(fn,"r") do file
        read(file, "hessian")
    end

    S_diag, U = eig(hessian)
    big_eig_vals = find(x -> x>1e-6, S_diag)
    rank = length(big_eig_vals)
    n = length(params)
    S_inv = zeros(n,n)
    for i = (n-rank+1):n
        S_inv[i,i] = 1/S_diag[i]
    end
    hessian_inv = U*sqrt(S_inv) #inverse of Hessian
    dist = DSGE.DegenerateMvNormal(params, hessian_inv)
 
    if :draw_type==mode 
        if parallel
            draws = pmap(x -> draw_from_approximated_post(m,data,dist), 1:n_part)
        else
            draws = [draw_from_approximated_post(m,data,dist) for j = 1:n_part]
        end
    else
        if parallel
            draws = pmap(x -> DSGE.draw_from_prior(m,data), 1:n_part)
        else
            draws = [DSGE.draw_from_prior(m,data)  for j = 1:n_part]
        end
end

    for i = 1:n_part
        prior_sim[i,:] = draws[i][1]
        logpost[i]     = draws[i][2]
        loglh[i]       = draws[i][3]
    end

    para_sim[1,:,:] = prior_sim # Draws from prior
    weight_sim[:,1] = 1/n_part
    zhat[1] = sum(weight_sim[:,1])

    i = 1
    para = para_sim[i, :, :]
    weight = repmat(weight_sim[:, i], 1, n_params)

    if m.testing
        open(rawpath(m,"estimate","draws.csv"),"a") do x
            writecsv(x,reshape(para_sim[i,:,:],(n_part,n_parameters(m))))
        end
    end

    μ = sum(para.*weight,1)
    σ = sqrt(sum((para - repmat(μ, n_part, 1)).^2 .*weight,1))


    if VERBOSITY[verbose] >= VERBOSITY[:low]
	println("--------------------------")
        println("Iteration = $(i) / $(n_Φ)")
	println("--------------------------")
        println("phi = $(tempering_schedule[i])")
	println("--------------------------")
        println("c = $(c)")
	println("--------------------------")
        if VERBOSITY[verbose] >= VERBOSITY[:high]
            for n=1:n_params
                println("$(m.parameters[n].key) = $(round(μ[n],5)), $(round(σ[n],5))")
	    end
        end
    end


    #RECURSION

    if VERBOSITY[verbose] >= VERBOSITY[:low]
	println("\n\n SMC Recursion starts \n\n")
    end

    total_sampling_time = 0.0

    for i=2:n_Φ
        tic()
	#------------------------------------
	# (a) Correction
	#------------------------------------
	# Incremental weights
	incweight = exp((tempering_schedule[i] - tempering_schedule[i-1]) * loglh)

       	# Update weights
	weight_sim[:,i] = weight_sim[:,i-1].*incweight
        zhat[i] = sum(weight_sim[:,i])

	# Normalize weights
	weight_sim[:, i] = weight_sim[:, i]/zhat[i]
	#------------------------------------
	# (b) Selection
	#------------------------------------

	ESS = 1/sum(weight_sim[:,i].^2)
        @assert !isnan(ESS) "no particles have non-zero weight"

	if (ESS < n_part/2)
	    id = resampler(weight_sim[:,i])
            # sample(1:n_part, weights(weight_sim[:,i]), n_part, replace = true)
            # systematic_resampling(m, weight_sim[:,i])
            para_sim[i-1, :, :] = para_sim[i-1, id, :]
	    loglh = loglh[id]
	    logpost = logpost[id]
	    weight_sim[:,i] = 1/n_part
	    n_resamp += 1
	    rsmp_sim[i] = 1
	end

	#------------------------------------
	# (c) Mutation
	#------------------------------------

        c = c*(0.95 + 0.10*exp(16*(accept - target))/(1 + exp(16*(accept - target))))
        m <= Setting(:step_size_smc, c)

        para = para_sim[i-1, :, free_para_inds]
        weight = repmat(weight_sim[:,i], 1, n_params)[:,free_para_inds]
        μ = sum(para.*weight,1)
	z =  (para - repmat(μ, n_part, 1))
        R_temp = (z.*weight)'*z
        R = nearest_spd(R_temp)

	temp_accept = zeros(n_part, 1) # Initialize acceptance indicator

        if parallel
            out = pmap(j -> mutation_RWMH(m, data, vec(para[j,:]'), loglh[j], logpost[j], i, R, tempering_schedule), 1:n_part)
            out = [mutation_RWMH(m, data, vec(para[j,:]'), loglh[j], logpost[j], i, R, tempering_schedule)  for j = 1:n_part]
        end

        para_sim_temp = Array{Float64}[out[i][1] for i=1:length(out)]
        fillmat = zeros(length(para_sim_temp), length(para_sim_temp[1]))
        for j = 1:length(para_sim_temp)
            fillmat[j,:] = para_sim_temp[j]
        end
        para_sim[i,:,:] = fillmat
        loglh = Float64[out[i][2] for i=1:length(out)]
        logpost = Float64[out[i][3] for i=1:length(out)]
        temp_accept = Int[out[i][4] for i=1:length(out)]

        accept = mean(temp_accept) # update average acceptance rate
        c_sim[i,:]    = c # scale parameter
        ESS_sim[i,:]  = ESS # ESS
        accept_sim[i,:] = accept # average acceptance rate

        para = para_sim[i, :, :]
        weight = repmat(weight_sim[:, i], 1, n_params)

        μ  = sum(para.*weight,1)
        σ = sqrt(sum((para - repmat(μ, n_part, 1)).^2 .*weight,1))

        # timekeeping
        block_time = toq()
        total_sampling_time += block_time
        total_sampling_time_minutes = total_sampling_time/60
        expected_time_remaining_sec = (total_sampling_time/i)*(n_Φ - i)
        expected_time_remaining_minutes = expected_time_remaining_sec/60

        if VERBOSITY[verbose] >= VERBOSITY[:low]
	    println("--------------------------")
            println("Iteration = $(i) / $(n_Φ)")
            println("time elapsed: $(round(total_sampling_time_minutes,4)) minutes")
            println("estimated time remaining: $(round(expected_time_remaining_minutes,4)) minutes")
	    println("--------------------------")
            println("phi = $(tempering_schedule[i])")
	    println("--------------------------")
            println("c = $(c)")
            println("accept = $(accept)")
            println("ESS = $(ESS)   ($(n_resamp) total resamples.)")
	    println("--------------------------")
            if VERBOSITY[verbose] >= VERBOSITY[:high]
                for n=1:n_params
                    println("$(m.parameters[n].key) = $(round(μ[n],5)), $(round(σ[n],5))")
	        end
            end
        end

    end
    if m.testing
        open(rawpath(m,"estimate,"wtsim.csv"),"a") do x
            writecsv(x,weight_sim)
        end
    end
    #------------------------------------
    # Saving parameter draws
    #------------------------------------
    if !m.testing

        save_file = h5open(rawpath(m,"estimate","smcsave.h5"),"w")
        n_saved_obs = n_part
        n_chunks = max(floor(n_saved_obs/5000),1)
        para_save = d_create(save_file, "smcparams", datatype(Float32),
                             dataspace(n_saved_obs,n_params),
                             "chunk", (n_chunks,n_params))

        print("\n saving output \n ")

        para_save[:,:] = map(Float32, para)

        close(save_file)
    end

end

function smc(m::AbstractModel, data::DataFrame; verbose::Symbol=:low)
    data_mat = df_to_matrix(m, data)
    return smc(m, data_mat, verbose=verbose)
end

function smc(m::AbstractModel; verbose::Symbol=:low)
    data = load_data(m)
    data_mat = df_to_matrix(m, data)
    return smc(m, data_mat, verbose=verbose)
end


"""
```
draw_from_prior(m::AbstractModel, data::Matrix)
```

Draw from prior distribution and resample if certain conditions are not met.
Returns a tuple (prior_draw, logpost, loglh).

"""
function draw_from_prior(m::AbstractModel, data::Matrix)
    success = false
    n_params = n_parameters(m)
    T = typeof(m.parameters[1].value)
    prior_draw = Array{T}(n_params)
    logpost = 0
    loglh = 0

    while !success
        for j in 1:n_params
            prior_draw[j] = if !m.parameters[j].fixed
                prior = rand(m.parameters[j].prior.value)

                # Resample until all prior draws are within the value bounds
                while !(m.parameters[j].valuebounds[1] < prior < m.parameters[j].valuebounds[2])
                    prior = rand(m.parameters[j].prior.value)
                end
                prior
            else
                m.parameters[j].value
            end
        end
        try
            logpost = posterior!(m, convert(Array{Float64,1},prior_draw), data)
            loglh = logpost - prior(m)
            success = true
        end
    end

    return (prior_draw, logpost, loglh)
end

