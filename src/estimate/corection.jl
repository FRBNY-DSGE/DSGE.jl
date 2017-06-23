function correction(s::Array{Float64}, ε::Array{Float64}, φ1:: Float64,φn::Float64, φn_1::Float64,w::Array{Float64},Φ::Array{Float64}, A::Array{Float64}, B::Array{Float64}, H::Array{Float64}, yt::Array{Float64,1})
    #s is 8xnumParticles
    #ε is 3xnumParticles
    #φn is scalar
    #φt1 is scalar (phi n-1)
    #w is a vector of weights
    #B is psi, the measurement error function (here since linear, we just premultiply it by s
    count = 2
    numParticles= length(s)
    check_InEff=InEff(1.0, φ1, yt, perror, H) 
    s_fore = Φ*s + sqrtS2*ε 
    perror = repmat(yt-A, 1, M)-B*s
    density=Array{numParticles,1}
    φ_old = φ1
    
    while (check_InEff>rstar)
        count = count + 1
        InEff_func(φ) = Ineff(φ, φ_old, yt, perror,H)-rstar
        ###the if prod sign thing seems uncessary. why can't just update check_InEff at end of while loop with the new φ before starting back around again?
        
        fphi_interval = [InEff_func(phi_old) InEff_func(1.0)]
        #Does this just mean subtract rstar???????? What are you doing, friend??
        if prod(sign(fphi_interval))==-1
            s_fore = Φ*s + sqrtS2*ε 
            perror = repmat(yt-A, 1, M)-B*s

            for i=1:numParticles
                density[i]= (φn/φ_old)^(length(yt)/2)*exp((-1/2)*perror[:,i]'*(φn-φ_old)*inv(H)*perror[:,i])
            end
            
            weights = (density.*w)./(sum(density.*weights))
            # call resampling method
            # We need to just have this in a for-loop... otherwise space-inefficient 
            weights = ones(numParticles,1)
            lik[t] = lik[t]+log(mean(density.*weights))
            
            c = c*(0.95 + 0.10*exp(16*(acpt-trgt))/(1+exp(16*(acpt-trgt)))
            #This syntax may be incorrect - either for getting the mean or what eps_up translates to
            μ = mean(ε_up,2)
            cov_s = (1/N)*(ε_up-repmat(μ,1,N))*(ε_up - repmat(μ,1,N))'
            checkp = diagm(diag(chol(A)))
           #WE ARE FEELING GOOD ABOUT THE PROSPECTS OF CHECKP BEING UPPER TRIANGULAR
           # if checkp ~= 0
           #    R = diag(diag(cov_s))
           # else
           #    R = cov_s
        else 
            check_InEff = rstar
        end

end