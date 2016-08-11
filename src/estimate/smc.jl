#Tune type to be refactored at a later date since the Tune type approach may not be the most
#intuitive given some of the fields are computational whereas some are theory driven. 
#Hence it's not the best way of organizing
"""
```
Tune
```

Tune contains the hyperparameters for "tuning" the SMC algorithm itself, inputs for tuning the MH algorithm, moments of the parameters to be given as output for the algorithm, and the indices of the fixed and non-fixed parameters for the purpose of covariance matrix construction.

### Fields

- `npara::Int64`: The number of parameters contained in m.parameters.
- `npart::Int64`: The number of particles being propogated by the algorithm.
- `nphi::Int64`: The number of stages for the tempering schedule.
- `lam::Float64`: The 'bending coefficient' λ in Φ(n) = (n/N(Φ))^λ.
- `phi::Array{Float64}`: The tempering schedule - an array with elements corresponding to the n stages of Φ.

- `nblock::Int64`: The number of random blocks in the Random Block Metropolis Hastings.
- `c::Float64`: The scaling factor for the covariance, constructed from acpt and trgt in the mutation step.
- `acpt::Float64`: The average acceptance rate of new particles in the mutation step.
- `trgt::Float64`: The target acceptance rate of particles.

- `mu::Array{Float64}`: The calculated mean of the parameters to be output after each recursion
- `R::Array{Float64}`: The calculated covariance matrix of the parameters to be output after each recursion.

- `f_inds::Array{Int64}`: An array of the indices of the fixed parameters
- `nf_inds::Array{Int64}`: An array of the indices of the non-fixed parameters
"""
type Tune
    npara::Int64 
    npart::Int64
    nphi::Int64
    lam::Float64
	phi::Array{Float64}

    nblock::Int64
    c::Float64
    acpt::Float64
    trgt::Float64
    
    mu::Array{Float64}
    R::Array{Float64}
    
	f_inds::Array{Int64}
    nf_inds::Array{Int64}
end

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

tune = Tune(length(m.parameters),1000,100,4.,[],6,0.5,0.25,0.25,[], [], [], [])

#Creating the tempering schedule
tune.phi = ((collect(1:1:tune.nphi)-1)/(tune.nphi-1)).^tune.lam

#Matrices for storing

parasim = zeros(tune.nphi, tune.npart, tune.npara) #parameter draws
wtsim = zeros(tune.npart, tune.nphi) #weights
zhat = zeros(tune.nphi, 1) #normalization constant
nresamp = 0 #record # of iteration resampled

csim = zeros(tune.nphi,1) #scale parameter
ESSsim = zeros(tune.nphi,1) #ESS
acptsim = zeros(tune.nphi,1) #average acceptance rate
rsmpsim = zeros(tune.nphi,1) #1 if resampled

#--------------------------------------------------------------
#Initialize Algorithm: Draws from prior
#--------------------------------------------------------------

if VERBOSITY[verbose] >= VERBOSITY[:low]
	println("\n\n SMC starts ....  \n\n  ")
end

#Particle draws from the parameter's marginal priors
priorsim = zeros(tune.npart,tune.npara)

# Posterior values at prior draws
loglh = zeros(tune.npart, 1)
logpost = zeros(tune.npart, 1)

for i in 1:tune.npart
    priodraw = []
    #Parameter draws per particle
    for j in 1:length(m.parameters)
        if !m.parameters[j].fixed
            prio = rand(m.parameters[j].prior.value)
            #Ensure all prior draws are within the value bounds
            while true
                if m.parameters[j].valuebounds[1] < prio && prio < m.parameters[j].valuebounds[2]
                    append!(priodraw, [prio])
                    break
                else
                    prio = rand(m.parameters[j].prior.value)
                end
            end
        else
            append!(priodraw, [m.parameters[j].value])
        end
    end
    priorsim[i,:] = priodraw'
    out = posterior!(m, convert(Array{Float64,1},priodraw), data; phi_smc = tune.phi[1])
    logpost[i] = out[:post]
    loglh[i] = out[:like]
end
parasim[1,:,:] = priorsim #Draws from prior #Lay priorsim draws on top of parasim box matrix which is 100x1000x13
wtsim[:,1] = 1/tune.npart #Initial weights are all equal, 1000x1
zhat[1] = round(sum(wtsim[:,1]),14) # zhat is 100x1 and its first entry is the sum of the first column of wtsim, the weights matrix

#RECURSION

tic()
totaltime = 0 #Probably let's take this out

if VERBOSITY[verbose] >= VERBOSITY[:low]
	println("\n\n SMC Recursion starts \n\n")
end

for i=2:tune.nphi
	#------------------------------------
	# (a) Correction
	#------------------------------------
	# Incremental weights
	incwt = exp((tune.phi[i] - tune.phi[i-1]) * loglh)

	# Update weights
	wtsim[:,i] = wtsim[:,i-1].*incwt #fill in other columns of wtsim
	
    zhat[i] = sum(wtsim[:,i]) # Fill in other entries of zhat with weights of each respective column of wtsim

	#Normalize weights
	wtsim[:, i] = wtsim[:, i]/zhat[i]
	#------------------------------------
	# (b) Selection 
	#------------------------------------

	ESS = 1/sum(wtsim[:,i].^2)
	if (ESS < tune.npart/2)
		id = systematic_resampling(wtsim[:,i])
        parasim[i-1, :, :] = parasim[i-1, id, :]
		loglh = loglh[id]
		logpost = logpost[id]
		wtsim[:,i] = 1/tune.npart
		nresamp += 1
		rsmpsim[i] = 1
	end

	#------------------------------------
	# (c) Mutation
	#------------------------------------

    tune.c = tune.c*(0.95 + 0.10*exp(16*(tune.acpt - tune.trgt))/(1 + exp(16*(tune.acpt - tune.trgt))))
   
    #For handling fixed parameters 
    if i ==2
        for k in 1:length(m.parameters)
            if m.parameters[k].fixed == true
                append!(tune.f_inds, [k])
            else
                append!(tune.nf_inds, [k])
            end
        end
    end

    para = squeeze(parasim[i-1, :, :][:,:,tune.nf_inds],1)
    wght = repmat(wtsim[:,i], 1, tune.npara)[:,tune.nf_inds]
    tune.mu = sum(para.*wght,1)
	z =  (para - repmat(tune.mu, tune.npart, 1))
    R_temp = (z.*wght)'*z
    tune.R = (R_temp + R_temp')/2

	#Particle mutation (RWMH 2)
	temp_acpt = zeros(tune.npart, 1) # Initialize acceptance indicator

    #Blocking RWMH
    #Generate random draws for each element in m.parameters
    #para_shuffled = the randomly shuffled vector of parameters (could be refactored, but
    #each code line follows the algorithm for random blocking suggested in the Schorfheide book)
    #block_vecs is a vector of tune.nblocks blocks, each block being a portion of the randomly
    #shuffled vector 

	#=
	"""
	`block(para::Array{Float64})`

	Following the procedure specified in Herbst & Schorfheide's "Bayesian Estimation of DSGE Models" to randomly shuffle and divide parameter vectors into blocks of similar/equal length given the input parameters vector from the model and the number of blocks from the tune object.

	### Argument:
	- `para`: A vector of parameters contained in m.parameters

	### Output:
	- `block_vecs`: An array of tune.nblocks elements (also an array), with each element being a "block" of the input vector para. 
	"""
	=#
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

    #for j = 1:tune.npart
        #block_para = block(para[j,:]')
        #calculate loglh per block
        #calculate logpost per block
        #for b = 1:tune.nblocks
        #ind_para, ind_loglh, ind_post, ind_acpt = mutation_RWMH(vec(block_para[b]), block_lh[b], block_post[b], tune, i, )
        #How to calculate the posteriors within the mutation step given data and the model are needed to calculate it. 
        #Would the data/model be modified to be blocked at the same divisions as the parameters?

    out = @sync @parallel (hcat) for j = 1:tune.npart 
        mutation_RWMH(vec(para[j,:]'), loglh[j], logpost[j], tune, i, data,m)
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
      
	tune.acpt = mean(temp_acpt) #update average acceptance rate	

    # store
    csim[i,:]    = tune.c; # scale parameter
    ESSsim[i,:]  = ESS; # ESS
    acptsim[i,:] = tune.acpt; # average acceptance rate

    # print some information (is this a useless conditional? i%1 will always == 0)
    if i % 1 == 0 
        
        para = squeeze(parasim[i, :, :],1);
        wght = repmat(wtsim[:, i], 1, length(m.parameters));

        mu  = sum(para.*wght,1);
        sig = sum((para - repmat(mu, tune.npart, 1)).^2 .*wght,1);
        sig = (sqrt(sig));

        # time calculation
	
	if VERBOSITY[verbose] == VERBOSITY[:high]
		println("--------------------------")
        println("Iteration = $(i) / $(tune.nphi)")
		println("--------------------------")
        println("phi = $(tune.phi[i])")
		println("--------------------------")
        println("c = $(tune.c)")
        println("acpt = $(tune.acpt)")
        println("ESS = $(ESS)   ($(nresamp) total resamples.)")
		println("--------------------------")
        for n=1:length(m.parameters)
            println("$(m.parameters[n].key) = $(mu[n]), $(sig[n])")
		end
    end
end
end
end
