#Using Debug

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
	c::AbstractFloat # of random blocks
	acpt::AbstractFloat # initial scale cov
	trgt::AbstractFloat #initial acpt rate
	alp::AbstractFloat #target acpt rate
	phi::Array{Float64}	#Mixture weight for mixture proposal

	R::Array # Covariance Matrix
	Rdiag::Array # Empty matrix but for the diagonal of R
	Rchol::Array # Lower triangular Cholesky Decomposition of R
	Rchol2::Array # Square root of Rdiag
end

#Instantiating a tune type (analogous to the struct in Matlab)
#The tempering schedule is created as the last argument in the constructor TuneType()
tune = Tune(length(m.parameters),1000,100,2,0.5,0.25,0.25,0.9,((collect(1:1:100)-1)/(100-1)).^2, [], [], [], [])

function smc(m::AbstractModel, data::Matrix, tune::Tune)
#--------------------------------------------------------------
#Set Parameters of Algorithm
#--------------------------------------------------------------

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
priorsim = zeros(tune.npart,tune.npara)
for i in 1:tune.npart
    priodraw = []
    for j in 1:length(m.parameters)
        #try catch for if a parameter is fixed
        try
            append!(priodraw, [rand(m.parameters[j].prior.value)])
        catch
            append!(priodraw, [m.parameters[j].value])
        end
    end
    priorsim[i,:] = priodraw'
end

parasim[1,:,:] = priorsim #Draws from prior #Lay priorsim draws on top of parasim box matrix which is 100x1000x13
wtsim[:,1] = 1/tune.npart #Initial weights are all equal, 1000x1
zhat[1] = sum(wtsim[:,1]) # zhat is 100x1 and its first entry is the sum of the first column of wtsim, the weights matrix

# Posterior values at prior draws
loglh = zeros(tune.npart, 1)
logpost = zeros(tune.npart, 1)

for i=1:1:tune.npart
	p0 = priorsim[i,:]';
	#logpost[i], loglh[i] = objfcn_dsge(p0, tune.phi[1], prio, bounds, data)
    logpost[i] = posterior(m, data; phi_smc = tune.phi[1])[:post]
    loglh[i] = likelihood(m,data)[1]
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
        parasim[i-1, :, :] = squeeze(parasim[i-1, id, :])
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
	@bp
    para = squeeze(parasim[i-1, :, :])
    wght = repmat(wtsim[:,i], 1, tune.npara)
	tune.mu = sum(para.*wght)
	z =  (para - repmat(tune.mu, tune.npart, 1))
	tune.R = (z.*wgth)'*z
	tune.Rdiag = diag(diag(tune.R))
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

	tune.acpt = mean(temp_acpt) #update average acceptance rat	


    # store
    csim[i,:]    = tune.c; # scale parameter
    ESSsim[i,:]  = ESS; # ESS
    acptsim[i,:] = tune.acpt; # average acceptance rate
    
    # print some information
    if mod[i, 1] == 0
        
        para = squeeze(parasim[i, :, :]);
        wght = repmat(wtsim[:, i], 1, 13);

mu  = sum(para.*wght);
sig = sum((para - repmat[mu, tune.npart, 1]).^2 .*wght);
sig = (sqrt(sig));

        # time calculation
		toc()

        #= To fix - string formatting

		println("--------------------------")
		println("Iteration = $i / tune.nphi")
		println("--------------------------")
		println("phi = tune.phi[i]")
		println("--------------------------")
		println("c = " tune.c)
		println("acpt = " tune.acpt)
		println("ESS = " tune.phi)
		println("--------------------------")
		param_names = ["tau", "kappa", "psil", "phi2", "rA", "piA", "gammaQ", "rho_R", "rho_G", "rho_z", "sigma_R", "sigma_g", "sigma_z"]
		for n=1:13
		    println(param_names[n]" = ", mu[n], sig[n])
		end
        =#
    end
end
end
