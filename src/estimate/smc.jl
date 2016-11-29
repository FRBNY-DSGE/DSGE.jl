#Blocking RWMH
#Generate random draws for each element in m.parameters
#para_shuffled = the randomly shuffled vector of parameters (could be refactored, but
#each code line follows the algorithm for random blocking suggested in the Schorfheide book)
#block_vecs is a vector of tune.nblocks blocks, each block being a portion of the randomly
#shuffled vector 
	
"""
`block(para::Array{Float64})`

Following the procedure specified in Herbst & Schorfheide's "Bayesian Estimation of DSGE Models" to randomly shuffle and divide parameter vectors into blocks of similar/equal length given the input parameters vector from the model and the number of blocks from the tune object.

### Argument:
- `para`: A vector of parameters contained in m.parameters

### Output:
- `block_vecs`: An array of tune.nblocks elements (also an array), with each element being a "block" of the input vector para. 
	"""
function block(para::Array{Float64})
    rand_order = [(n,rand()) for n in 1:length(para)]
    rand_tups = sort(rand_order, by = x -> x[2])
    rand_inds = [n[1] for n in rand_tups]
    para_shuffled = para[rand_inds]
    
    if length(para)%tune.nblocks == 0
        block_size = length(para)/tune.nblocks
        append!(block_size_iter, [block_size for _ in 1:tune.nblocks])   
    else
        block_size = round(length(para)/tune.nblocks)+1
        append!(block_size_iter, [block_size-1 for _ in 1:block_size*tune.nblocks-length(para)])
        append!(block_size_iter, [block_size for _ in 1:tune.nblocks-(block_size*tune.nblocks-length(para))])
        
    end
    
    block_size_iter = []
    start_block = 1
    end_block = block_size_iter[1]
    
    #off = the off-indices of when the blocks change 
    off = [0,cumsum(block_size_iter)]
    block_vecs = [para_shuffled[off[i]+1 : off[i]+block_size_iter[i]]
                  for i=1:tune.nblocks]
    return block_vecs
end

#for j = 1:n_part
#block_para = block(para[j,:]')
#calculate loglh per block
#calculate logpost per block
#for b = 1:tune.nblocks
#ind_para, ind_loglh, ind_post, ind_acpt = mutation_RWMH(vec(block_para[b]), block_lh[b], block_post[b], tune, i, )
#How to calculate the posteriors within the mutation step given data and the model are needed to calculate it. 
#Would the data/model be modified to be blocked at the same divisions as the parameters?


"""
```
smc(m,data)
```

### Inputs

- `m`: A model object, from which its parameters values, prior dists, and bounds will be referenced
- `data`: A matrix containing time series of the observables to be used in the calculation of the posterior/likelihood

### Optional Inputs
- `verbose`: The desired frequency of function progress messages printed to standard out.
	- `:none`: No status updates will be reported.
	- `:low`: Status updates for SMC initialization and recursion will be included.
	    - `:high`: Status updates for every iteration of SMC is output, which includes the mu and sigma of each individual parameter draw after each iteration, as well as the calculated acceptance rate, the ESS, and the number of times resampled.

### Outputs

- `mu`: The mean of a particular parameter as calculated by taking the weighted sum of that parameter from the particle draws made in that particular stage Nφ. 
- `sig`: The standard deviation of a particular parameter. 

### Overview

Sequential Monte Carlo can be used in lieu of Random Walk Metropolis Hastings to generate parameter samples from high-dimensional parameter spaces using sequentially constructed proposal densities.  

SMC is broken up into three main steps:

- `Correction`: Reweight the particles from stage n-1 by defining "incremental weights", incwt, which gradually incorporate the likelihood function p(Y|θ(i,n-1)) into the particle weights. 
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
    c = get_setting(m, :c)
    acpt = get_setting(m, :acpt)
    trgt = get_setting(m, :trgt)

    parallel = get_setting(m, :use_parallel_workers)

    #Creating the tempering schedule
    λ = get_setting(m, :λ)
    tempering_schedule = ((collect(1:1:n_Φ)-1)/(n_Φ-1)).^λ

    #Matrices for storing

    parasim = zeros(n_Φ, n_part, n_params) #parameter draws
    wtsim = zeros(n_part, n_Φ) #weights
    zhat = zeros(n_Φ, 1) #normalization constant
    nresamp = 0 #record # of iteration resampled

    csim = zeros(n_Φ,1) #scale parameter
    ESSsim = zeros(n_Φ,1) #ESS
    acptsim = zeros(n_Φ,1) #average acceptance rate
    rsmpsim = zeros(n_Φ,1) #1 if resampled

    #--------------------------------------------------------------
    #Initialize Algorithm: Draws from prior
    #--------------------------------------------------------------

    if VERBOSITY[verbose] >= VERBOSITY[:low]
	println("\n\n SMC starts ....  \n\n  ")
    end

    #Particle draws from the parameter's marginal priors
    priorsim = zeros(n_part,n_params)

    # Posterior values at prior draws
    loglh = zeros(n_part, 1)
    logpost = zeros(n_part, 1)

    # fix seed if testing
    if m.testing
        println("seed set")
        srand(42)
    end

    for i in 1:n_part
        T = typeof(m.parameters[1].value)
        priodraw = Array{T}(n_params)

        #Parameter draws per particle
        for j in 1:n_params

            priodraw[j] = if !m.parameters[j].fixed            
                prio = rand(m.parameters[j].prior.value)
                
                # Resample until all prior draws are within the value bounds
                while !(m.parameters[j].valuebounds[1] < prio < m.parameters[j].valuebounds[2])
                    prio = rand(m.parameters[j].prior.value)
                end

                prio
            else
                m.parameters[j].value
            end
        end
        priorsim[i,:] = priodraw'
        out = posterior!(m, convert(Array{Float64,1},priodraw), data; phi_smc = tempering_schedule[1])
        logpost[i] = out[:post]
        loglh[i] = out[:like]
    end
    parasim[1,:,:] = priorsim #Draws from prior #Lay priorsim draws on top of parasim box matrix which is 100x1000x13
    wtsim[:,1] = 1/n_part #Initial weights are all equal, 1000x1
    zhat[1] = round(sum(wtsim[:,1]),14) # zhat is 100x1 and its first entry is the sum of the first column of wtsim, the weights matrix
        
    # print info for prior draws
    i = 1
    para = squeeze(parasim[i, :, :],1);
    wght = repmat(wtsim[:, i], 1, n_params);
    
    if m.testing
        open("draws.csv","a") do x
            writecsv(x,reshape(parasim[i,:,:],(n_part,n_parameters(m))))
        end 
    end
     
    mu  = sum(para.*wght,1);
    sig = sum((para - repmat(mu, n_part, 1)).^2 .*wght,1);
    sig = (sqrt(sig));
    
    # time calculation
    
    if VERBOSITY[verbose] >= VERBOSITY[:low]
	println("--------------------------")
        println("Iteration = $(i) / $(n_Φ)")
	println("--------------------------")
        println("phi = $(tempering_schedule[i])")
	println("--------------------------")
        println("c = $(c)")
	println("--------------------------")
        for n=1:n_params
            println("$(m.parameters[n].key) = $(mu[n]), $(sig[n])")
	end
    end


    #RECURSION

    tic()
    totaltime = 0 #Probably let's take this out

    if VERBOSITY[verbose] >= VERBOSITY[:low]
	println("\n\n SMC Recursion starts \n\n")
    end

    for i=2:n_Φ
	#------------------------------------
	# (a) Correction
	#------------------------------------
	# Incremental weights
	incwt = exp((tempering_schedule[i] - tempering_schedule[i-1]) * loglh)

	# Update weights
	wtsim[:,i] = wtsim[:,i-1].*incwt #fill in other columns of wtsim
	
        zhat[i] = sum(wtsim[:,i]) # Fill in other entries of zhat with weights of each respective column of wtsim

	#Normalize weights
	wtsim[:, i] = wtsim[:, i]/zhat[i]
	#------------------------------------
	# (b) Selection 
	#------------------------------------

	ESS = 1/sum(wtsim[:,i].^2)
	if (ESS < n_part/2)
	    id = systematic_resampling(m, wtsim[:,i])
            parasim[i-1, :, :] = parasim[i-1, id, :]
	    loglh = loglh[id]
	    logpost = logpost[id]
	    wtsim[:,i] = 1/n_part
	    nresamp += 1
	    rsmpsim[i] = 1
	end

	#------------------------------------
	# (c) Mutation
	#------------------------------------

        c = c*(0.95 + 0.10*exp(16*(acpt - trgt))/(1 + exp(16*(acpt - trgt))))
        
        para = squeeze(parasim[i-1, :, :][:,:,free_para_inds],1)
        wght = repmat(wtsim[:,i], 1, n_params)[:,free_para_inds]
        mu = sum(para.*wght,1)
	z =  (para - repmat(mu, n_part, 1))
        R_temp = (z.*wght)'*z
        R = (R_temp + R_temp')/2

	#Particle mutation (RWMH 2)
	temp_acpt = zeros(n_part, 1) # Initialize acceptance indicator
        
        if parallel
            out = @sync @parallel (hcat) for j = 1:n_part 
                mutation_RWMH(m, data, vec(para[j,:]'), loglh[j], logpost[j], i, R, tempering_schedule)
            end
        else
            out = [mutation_RWMH(m, data, vec(para[j,:]'), loglh[j], logpost[j], i, R, tempering_schedule)  for j = 1:n_part]
        end

        parasimtemp = Array{Float64}[out[i][1] for i=1:length(out)]
        fillmat = zeros(length(parasimtemp), length(parasimtemp[1]))
        for j = 1:length(parasimtemp)
            fillmat[j,:] = parasimtemp[j]
        end
        parasim[i,:,:] = fillmat
        loglh = Float64[out[i][2] for i=1:length(out)]
        logpost = Float64[out[i][3] for i=1:length(out)]
        temp_acpt = Int[out[i][4] for i=1:length(out)]
        
        acpt = mean(temp_acpt) # update average acceptance rate
        csim[i,:]    = c; # scale parameter
        ESSsim[i,:]  = ESS; # ESS
        acptsim[i,:] = acpt; # average acceptance rate
        
        para = squeeze(parasim[i, :, :],1);
        wght = repmat(wtsim[:, i], 1, n_params);
        
        mu  = sum(para.*wght,1);
        sig = sum((para - repmat(mu, n_part, 1)).^2 .*wght,1);
        sig = (sqrt(sig));
        
        # time calculation        
        
        if VERBOSITY[verbose] >= VERBOSITY[:low]
	    println("--------------------------")
            println("Iteration = $(i) / $(n_Φ)")
	    println("--------------------------")
            println("phi = $(tempering_schedule[i])")
	    println("--------------------------")
            println("c = $(c)")
            println("acpt = $(acpt)")
            println("ESS = $(ESS)   ($(nresamp) total resamples.)")
	    println("--------------------------")
            if VERBOSITY[verbose] >= VERBOSITY[:high]
                for n=1:n_params
                    println("$(m.parameters[n].key) = $(mu[n]), $(sig[n])")
	        end
            end
        end
        
        
    end
    if m.testing
        open("wtsim.csv","a") do x
            writecsv(x,wtsim)
        end
    end
    #------------------------------------
    # Saving parameter draws
    #------------------------------------
    if !m.testing
        save_file = h5open(rawpath(m,"estimate","mhsave.h5"),"w")
        n_saved_obs = n_part
        n_chunks = max(floor(n_saved_obs/5000),1)
        para_save = d_create(save_file, "mhparams", datatype(Float32),
                             dataspace(n_saved_obs,n_params),
                             "chunk", (n_chunks,n_params))
        post_save = d_create(save_file, "mhposterior", datatype(Float32),
                             dataspace(n_saved_obs,1),
                             "chunk", (n_chunks,1))
        TTT_save  = d_create(save_file, "mhTTT", datatype(Float32),
                             dataspace(n_saved_obs,n_states_augmented(m),n_states_augmented(m)),
                             "chunk", (n_chunks,n_states_augmented(m),n_states_augmented(m)))
        RRR_save  = d_create(save_file, "mhRRR", datatype(Float32),
                             dataspace(n_saved_obs,n_states_augmented(m),n_shocks_exogenous(m)),
                             "chunk", (n_chunks,n_states_augmented(m),n_shocks_exogenous(m)))
        zend_save = d_create(save_file, "mhzend", datatype(Float32),
                             dataspace(n_saved_obs,n_states_augmented(m)),
                             "chunk", (n_chunks,n_states_augmented(m)))
        
        if parallel
            out = @sync @parallel (hcat) for j = 1:n_part 
                posterior!(m, vec(para[j,:]'), data, mh=true)
            end
        else
            out = [posterior!(m, vec(para[j,:]'), data, mh=true)  for j = 1:n_part]
        end
        
        for i = 1:n_part
            para_save[i,:]  = map(Float32,para[i,:])
            post_save[i,:]  = map(Float32,out[i][:post])
            TTT_save[i,:,:] = map(Float32,out[i][:mats][:TTT])
            RRR_save[i,:,:] = map(Float32,out[i][:mats][:RRR])
            zend_save[i,:]  = map(Float32,out[i][:mats][:zend]')
        end
        close(save_file)
    end
        
end
