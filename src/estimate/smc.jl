using Debug

#to be refactored at a later date, npara can be removed completely
#npart through phi could be written as optional arguments
#R's can be specified in the function itself

type Tune
	#General
	npara::Int64 # of parameters
	npart::Int64 # of particles
	nphi::Int64 # of stages
	lam::Int64 # bending coeff, lam = 1 means linear cooling schedule

	#Tuning for MH algorithms
	c::Float64 # of random blocks
	acpt::Float64 # initial scale cov
	trgt::Float64 #initial acpt rate
	alp::Float64 #target acpt rate
	phi::Array{Float64}	#Mixture weight for mixture proposal
    
    mu::Array{Float64} #Mean
    R::Array{Float64} # Covariance Matrix
    Rdiag::Array{Float64} # Empty matrix but for the diagonal of R
    Rchol::Array{Float64} # Lower triangular Cholesky Decomposition of R
    Rchol2::Array{Float64} # Square root of Rdiag
end



@debug function smc(m::AbstractModel, data::Matrix)
#--------------------------------------------------------------
#Set Parameters of Algorithm
#--------------------------------------------------------------


#Instantiating a tune type (analogous to the struct in Matlab)
#The tempering schedule is created as the last argument in the constructor TuneType()
tune = Tune(length(m.parameters),1000,100,2,0.5,0.25,0.25,0.9,((collect(1:1:100)-1)/(100-1)).^2,[], [], [], [], [])

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

println("\n\n SMC starts ....  \n\n  ")
#Draws from the prior
@bp
priorsim = zeros(tune.npart,tune.npara)
for i in 1:tune.npart
    priodraw = []
    for j in 1:length(m.parameters)
        #try catch for if a parameter is fixed
        if !m.parameters[j].fixed
            prio = rand(m.parameters[j].prior.value)
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
end
@bp
parasim[1,:,:] = priorsim #Draws from prior #Lay priorsim draws on top of parasim box matrix which is 100x1000x13
wtsim[:,1] = 1/tune.npart #Initial weights are all equal, 1000x1
zhat[1] = round(sum(wtsim[:,1]),14) # zhat is 100x1 and its first entry is the sum of the first column of wtsim, the weights matrix

# Posterior values at prior draws
loglh = zeros(tune.npart, 1)
logpost = zeros(tune.npart, 1)
for i=1:1:tune.npart
	p0 = priorsim[i,:]';
    #logpost[i], loglh[i] = objfcn_dsge(p0, tune.phi[1], prio, bounds, data)
    logpost[i] = posterior!(m, squeeze(p0,2), data; phi_smc = tune.phi[1])[:post]
    loglh[i] = posterior!(m, squeeze(p0,2), data; phi_smc = tune.phi[1])[:like]
end

#RECURSION

tic()
totaltime = 0 #Probably let's take this out

println("\n\n SMC Recursion starts \n\n")

for i=2:1:tune.nphi
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
    sampled = false
	if (ESS < tune.npart/2)
        sampled = true
		id, m = systematic_resampling(wtsim[:,i]')
        parasim[i-1, :, :] = squeeze(parasim[i-1, id, :],1)
		loglh = loglh[id]
		logpost = logpost[id]
		wtsim[:,i] = 1/tune.npart
		nresamp += 1
		rsmpsim[i] = 1
        if m.testing
            if !sampled
                h5open(workpath(m, "estimate","degen_dist.h5"),"w") do f
                    f["wtsim"] = wtsim[:,i]
                end
                sampled = true
            end
        end
	end

	#------------------------------------
	# (c) Mutation
	#------------------------------------
    para = squeeze(parasim[i-1, :, :],1)
    wght = repmat(wtsim[:,i], 1, tune.npara)
    tune.mu = sum(para.*wght,1)
	z =  (para - repmat(tune.mu, tune.npart, 1))
	tune.R = (z.*wght)'*z
    tune.Rdiag = diagm(diag(tune.R))
	tune.Rchol = chol(tune.R, Val{:L})
	tune.Rchol2 = sqrt(tune.Rdiag)

	#Particle mutation (RWMH 2)
	temp_acpt = zeros(tune.npart, 1) # Initialize acceptance indicator

	@parallel for j=1:tune.npart
		ind_para, ind_loglh, ind_post, ind_acpt = mutation_RMWH(para[j,:]', loglh[j], logpost[j], tune, i, f)
		parasim[i,j,:] = ind_para
		loglh[j] = ind_loglh
		logpost[j] = ind_post
		temp_acpt[j,1] = ind_acpt
	end

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
		toc()

		println("--------------------------")
        println("Iteration = $(i) / $(tune.nphi)")
		println("--------------------------")
        println("phi = $(tune.phi[i])")
		println("--------------------------")
        println("c = $(tune.c)")
        println("acpt = $(tune.acpt)")
        println("ESS = $(tune.phi)")
		println("--------------------------")
        for n=1:length(m.parameters)
            println("$(m.parameters[n].key) = ,$(mu[n]), $(sig[n])")
		end
    end
end
end
