"""
```
steadystate!(m::TwoAssetHANK)
```

Calculates the model's steady-state values. `steadystate!(m)` must be called whenever
the parameters of `m` are updated.
"""
function steadystate!(m::TwoAssetHANK)
    # Read in parameters
    aalpha = m[:aalpha].value
    ddelta = m[:ddelta].value
    ddeath = m[:ddeath].value
    rrho   = m[:rrho].value
    chi0   = m[:chi0].value
    chi1   = m[:chi1].value
    chi2   = m[:chi2].value
    a_lb   = m[:a_lb].value

    pam = m[:pam].value
    xxi = m[:xxi].value
    ggamma = m[:ggamma].value
    tau_I = m[:tau_I].value
    trans = m[:trans].value

    # Set liquid rates
    r_b_SS = m[:r_b_SS]
    r_b_borr_SS = m[:r_b_borr_SS]
    borrwedge_SS = m[:borrwedge_SS]

    K_liquid = get_setting(m, :K_liquid)
    aggregate_variables = get_setting(m, :aggregate_variables)
    distributional_variables = get_setting(m, :distributional_variables)
    distributional_variables_1 = get_setting(m, :distributional_variables_1)

    # Read in approximation parameters
    maxit_HJB  = get_setting(m, :maxit_HJB)
    crit_HJB   = get_setting(m, :crit_HJB)
    Delta      = get_setting(m, :Delta)

    maxit_HIS  = get_setting(m, :maxit_HIS)
    crit_HIS   = get_setting(m, :crit_HIS)
    start_HIS  = get_setting(m, :start_HIS)

    maxit_KFE  = get_setting(m, :maxit_KFE)
    crit_KFE   = get_setting(m, :crit_KFE)
    Delta_KFE  = get_setting(m, :Delta_KFE)

    maxit_KL   = get_setting(m, :maxit_KL)
    crit_KL    = get_setting(m, :crit_KL)
    relax_KL   = get_setting(m, :relax_KL)

    # Read in grids
    I       = get_setting(m, :I)
    J       = get_setting(m, :J)
    a_g     = get_setting(m, :a_g)
    b_g     = get_setting(m, :b_g)
    I_g     = get_setting(m, :I_g)
    J_g     = get_setting(m, :J_g)
    N       = get_setting(m, :N)
    a       = get_setting(m, :a)
    b       = get_setting(m, :b)
    y       = get_setting(m, :y)
    y_dist  = get_setting(m, :y_dist)
    y_mean  = get_setting(m, :y_mean)
    KL      = get_setting(m, :KL_0)

    interp_decision = get_setting(m, :interp_decision)

    # Set liquid rates
    r_b      = r_b_SS
    r_b_borr = r_b_borr_SS

    # Compute prices associated with initial guess of KL
    w	= (1 - aalpha) * (KL ^ aalpha)
    r_a	= aalpha * (KL ^ (aalpha - 1)) - ddelta

    a_grid, a_g_grid, b_grid, y_grid, y_g_grid, r_a_grid, r_b_grid, r_a_g_grid, r_b_g_grid, daf_grid, daf_g_grid, dab_grid, dab_g_grid, dab_g_tilde_grid, dbf_grid, dbf_g_grid, dbb_grid, dbb_g_grid, trans_grid, l_grid, l_g_grid, w_grid = set_grids(a, b, a_g, b_g, y, I, J, I_g, J_g, N, w, r_a, r_b, r_b_borr, trans)

    # Initial consumption and value function
    c_0	= (1-xxi) * w * y_grid .* l_grid + (r_a_grid[1,1,1] + ddeath*pam) .* a_grid +
        (r_b_borr + ddeath * pam) .* b_grid + trans_grid - tau_I * w * y_grid .* l_grid
    V_0	= 1/(rrho + ddeath) * u_fn(c_0, ggamma)

    # Initial distribution
    gg0 = zeros(I_g,J_g,N)
    gg0[b .== 0, 1, :] = y_dist
    gg0 = gg0 ./ sum(gg0)
    gg0 = gg0 ./ dab_g_tilde_grid 	# ensure integration to 1
    gg0 = reshape(gg0, I_g*J_g*N, 1)
    gg = gg0

    # RECANOW V_n =

    #----------------------------------------------------------------
    # Iterate on KL to find steady state
    #----------------------------------------------------------------
    for ii	= 1 : maxit_KL

	    # Derive aggregates given KL
	    w			= (1 - aalpha) * (KL ^ aalpha)
 	    r_a			= aalpha * (KL ^ (aalpha - 1)) - ddelta
	    r_a_grid    = repeat([r_a],I,J,N)
        r_a_g_grid	= repeat([r_a],I_g,J_g,N)
        w_g_grid	= repeat([w],I_g,J_g,N)

	    # Store current value function
	    Vn	= V_0

		# Preparations
		cF  = Array{Float64}(undef, I,J,N)
		sF  = Array{Float64}(undef, I,J,N)
		HcF = Array{Float64}(undef, I,J,N)

		cB  = Array{Float64}(undef, I,J,N)
		sB  = Array{Float64}(undef, I,J,N)
		HcB = Array{Float64}(undef, I,J,N)

		c0  = Array{Float64}(undef, I,J,N)
		s0  = Array{Float64}(undef, I,J,N)
		Hc0 = Array{Float64}(undef, I,J,N)

    	# Preparations
		dFB 	= Array{Float64}(undef, I,J,N)
		HdFB 	= Array{Float64}(undef, I,J,N)

		dBF		= Array{Float64}(undef, I,J,N)
		HdBF	= Array{Float64}(undef, I,J,N)

		dBB   	= Array{Float64}(undef, I,J,N)
		HdBB 	= Array{Float64}(undef, I,J,N)

	    #----------------------------------------------------------------
	    # Solve HJB
	    #----------------------------------------------------------------
        for nn    = 1 : maxit_HJB
            #-----
            # Compute derivatives w.r.t. illiquid assets a
            #-----
            # Preparations
            VaF        = zeros(I,J,N)
            VaB     = zeros(I,J,N)
            Vamin     = 0

            # Forward difference
            VaF[:,1:J-1,:]     = (Vn[:,2:J,:]-Vn[:,1:J-1,:]) ./ daf_grid[:,1:J-1,:]
            VaF[:,1:J-1,:]     = max.(VaF[:,1:J-1,:], Vamin)

            # Backward difference
            VaB[:,2:J,:]     = (Vn[:,2:J,:]-Vn[:,1:J-1,:]) ./ dab_grid[:,2:J,:]
            VaB[:,2:J,:]     = max.(VaB[:,2:J,:], Vamin)

            #------------------------------------------------------------
            # Compute derivatives w.r.t. liquid assets b
            #------------------------------------------------------------

            # Preparations
            VbF = zeros(I,J,N)
            VbB = zeros(I,J,N)
            Vbmin = 1e-8

            # Forward difference
            VbF[1:I-1,:,:] = (Vn[2:I,:,:]-Vn[1:I-1,:,:]) ./ dbf_grid[1:I-1,:,:]
            VbF[1:I-1,:,:] = max.(VbF[1:I-1,:,:],Vbmin)

            # Backward difference
            VbB[2:I,:,:] = (Vn[2:I,:,:]-Vn[1:I-1,:,:]) ./ dbb_grid[2:I,:,:]
            VbB[2:I,:,:] = max.(VbB[2:I,:,:],Vbmin)

            #------------------------------------------------------------
            # Consumption decisions
            #------------------------------------------------------------
            # Decisions conditional on a particular direction of derivative
            cF[1:I-1,:,:]     = VbF[1:I-1,:,:] .^ (-1/ggamma)
            cF[I,:,:]         = zeros(1,J,N)
            sF[1:I-1,:,:]     = ((1-xxi)-tau_I) * w * l_grid[1:I-1,:,:] .* y_grid[1:I-1,:,:] .+
                b_grid[1:I-1,:,:] .* (r_b_grid[1:I-1,:,:] .+ ddeath*pam) .+
                trans_grid[1:I-1,:,:] .- cF[1:I-1,:,:]
            sF[I,:,:]         = zeros(1,J,N)
            HcF[1:I-1,:,:]    = u_fn(cF[1:I-1,:,:],ggamma) .+ VbF[1:I-1,:,:] .* sF[1:I-1,:,:]
            HcF[I,:,:]         = -1e12*ones(1,J,N)
            validF            = (sF .> 0)

            cB[2:I,:,:]     = VbB[2:I,:,:].^(-1/ggamma)
            cB[1,:,:]         = ((1-xxi)-tau_I) * w * l_grid[1,:,:] .* y_grid[1,:,:] .+
                b_grid[1,:,:] .* (r_b_grid[1,:,:] .+ ddeath*pam) .+ trans_grid[1,:,:]
            sB[2:I,:,:]     = ((1-xxi)-tau_I) * w * l_grid[2:I,:,:] .* y_grid[2:I,:,:] .+
                b_grid[2:I,:,:] .* (r_b_grid[2:I,:,:] .+ ddeath*pam) .+ trans_grid[2:I,:,:] .-
                cB[2:I,:,:]
            sB[1,:,:]         = zeros(1,J,N)
            HcB[:,:,:]         = u_fn(cB, ggamma) .+ VbB .* sB
            validB            = (sB .< 0)

            c0[:,:,:]         = ((1-xxi)-tau_I) * w * l_grid[:,:,:] .* y_grid[:,:,:] .+
                b_grid[:,:,:] .* (r_b_grid[:,:,:] .+ ddeath*pam) .+ trans_grid[:,:,:]
            s0[:,:,:]         = zeros(I,J,N)
            Hc0[:,:,:]         = u_fn(c0,ggamma)

            # Which direction to use
            IcF             = validF .* max.(.!validB,(HcF.>=HcB)) .* (HcF.>=Hc0)
            IcB             = validB .* max.(.!validF,(HcB.>=HcF)) .* (HcB.>=Hc0)
            Ic0             = 1 .- IcF .- IcB

            # Optimal consumption and liquid savings
            c        = IcF .* cF + IcB .* cB + Ic0 .* c0
            s        = IcF .* sF + IcB .* sB + Ic0 .* s0
            u        = u_fn(c,ggamma)

            #------------------------------------------------------------
            # Deposit decision
            #------------------------------------------------------------
            # Decisions conditional on a particular direction of derivative
            dFB[2:I,1:J-1,:]     = opt_deposits(VaF[2:I,1:J-1,:], VbB[2:I,1:J-1,:],
                                                a_grid[2:I,1:J-1,:], chi0, chi1, chi2, a_lb)
            dFB[:,J,:]             = zeros(I,1,N)
            dFB[1,1:J-1,:]         = zeros(1,J-1,N)
            HdFB[2:I,1:J-1,:]     = VaF[2:I,1:J-1,:] .* dFB[2:I,1:J-1,:] - VbB[2:I,1:J-1,:] .*
                (dFB[2:I,1:J-1,:] +
                 adj_cost_fn(dFB[2:I,1:J-1,:], a_grid[2:I,1:J-1,:], chi0, chi1, chi2, a_lb))
            HdFB[:,J,:]         = -1.0e12 * ones(I,1,N)
            HdFB[1,1:J-1,:]     = -1.0e12 * ones(1,J-1,N)
            validFB             = (dFB .> 0) .* (HdFB .> 0)

            dBF[1:I-1,2:J,:]     = opt_deposits(VaB[1:I-1,2:J,:],VbF[1:I-1,2:J,:],
                                                a_grid[1:I-1,2:J,:], chi0, chi1, chi2, a_lb)
            dBF[:,1,:]    = zeros(I,1,N)
            dBF[I,2:J,:]  = zeros(1,J-1,N)
            HdBF[1:I-1,2:J,:]  = VaB[1:I-1,2:J,:] .* dBF[1:I-1,2:J,:] - VbF[1:I-1,2:J,:] .*
                (dBF[1:I-1,2:J,:] +
                 adj_cost_fn(dBF[1:I-1,2:J,:], a_grid[1:I-1,2:J,:], chi0, chi1, chi2, a_lb))
            HdBF[:,1,:]   = -1.0e12 * ones(I,1,N)
            HdBF[I,2:J,:]         = -1.0e12 * ones(1,J-1,N)
            validBF             = (dBF .<= -adj_cost_fn(dBF, a_grid, chi0, chi1, chi2, a_lb)) .*
                (HdBF .> 0)

            VbB[1,2:J,:]         = u_fn(cB[1,2:J,:],ggamma)
            dBB[:,2:J,:]         = opt_deposits(VaB[:,2:J,:],VbB[:,2:J,:],a_grid[:,2:J,:],
                                                chi0, chi1, chi2, a_lb)
            dBB[:,1,:]     = zeros(I,1,N)
            HdBB[:,2:J,:]  = VaB[:,2:J,:] .* dBB[:,2:J,:] - VbB[:,2:J,:] .*
                (dBB[:,2:J,:] + adj_cost_fn(dBB[:,2:J,:],a_grid[:,2:J,:], chi0, chi1, chi2, a_lb))
            HdBB[:,1,:]    = -1.0e12 * ones(I,1,N)
            validBB        = (dBB .> -adj_cost_fn(dBB, a_grid, chi0, chi1, chi2, a_lb)) .*
                (dBB .<= 0) .* (HdBB .> 0)

            # Which direction to use
            IcFB = validFB .* max.(.!validBF,(HdFB .>= HdBF)) .* max.(.!validBB,(HdFB .>= HdBB))
            IcBF = max.(.!validFB,(HdBF .>= HdFB)) .* validBF .* max.(.!validBB,(HdBF .>= HdBB))
            IcBB = max.(.!validFB,(HdBB .>= HdFB)) .* max.(.!validBF,(HdBB .>= HdBF)) .* validBB
            Ic00 = (.!validFB) .* (.!validBF) .* (.!validBB)

            # Optimal deposits
            d = IcFB .* dFB + IcBF .* dBF + IcBB .* dBB + Ic00 .* zeros(I,J,N)

            #------------------------------------------------------------
            # Construct "A" matrix multiplying value function
            #------------------------------------------------------------

            # interpolate
            d_g = reshape(interp_decision*d[:],I_g,J_g,N)
            s_g = reshape(interp_decision*s[:],I_g,J_g,N)
            c_g = reshape(interp_decision*c[:],I_g,J_g,N)

            aa, bb, aau, bbu = transition(I_g, J_g, N, I, J, ddeath, pam, xxi, w, chi0,
                                          chi1, chi2, a_lb, l_grid, l_g_grid, y_grid, y_g_grid,
                                          d, dab_grid, daf_grid, d_g,a_grid, a_g_grid, s, s_g,
                                          r_a_grid, r_a_g_grid)
            cc  = kron(lambda,my_speye(I*J))
            ccu = kron(lambda,my_speye(I_g*J_g))

            A = aa + bb

        #------------------------------------------------------------
        # Update value function
        #------------------------------------------------------------

        global Vn1 			= Array{Float64}(undef, I,J,N)
        global Bik_all 	    = Vector{Array{Float64}}(undef, N)#cell(N,1)
        global uk_stacked, Vk_stacked = Vector{Float64}(undef, I*J), Vector{Float64}(undef, I*J)
        for kk = 1:N
        	Ak 				= A[1+(kk-1)*(I*J):kk*(I*J),1+(kk-1)*(I*J):kk*(I*J)]
        	Bk 				= (1 + Delta*(rrho + ddeath) - Delta*lambda[kk,kk])*my_speye(I*J) - Delta*Ak
        	global Bik_all[kk] 	= inv(Matrix{Float64}(Bk))
        	global uk_stacked 		= reshape(u[:,:,kk],I*J,1)
            global Vk_stacked 		= reshape(Vn[:,:,kk],I*J,1)
        	indx_k 			= 1:N .!= kk    #~ismember(1:N,kk)
        	Vkp_stacked 	= sum(repeat(lambda[kk,indx_k]',I*J,1) .*
                                  reshape(Vn[:,:,indx_k],I*J,N-1), dims=2)
        	qk 				= Delta*uk_stacked + Vk_stacked + Delta*Vkp_stacked
        	Vn1k_stacked 	= Bik_all[kk]*qk
        	global Vn1[:,:,kk] 	= reshape(Vn1k_stacked,I,J,1)
        end

        # Howard improvement step
        if nn >= start_HIS
        	for jj = 1:maxit_HIS
        		Vn2 = Array{Float64}(undef, I, J, N)
        		for kk = 1:N
        			uk_stacked 	= reshape(u[:,:,kk],I*J,1)
		        	Vk_stacked 	= reshape(Vn1[:,:,kk],I*J,1)
        			indx_k 		= 1:N .!= kk
		        	Vkp_stacked	= sum(repeat(lambda[kk,indx_k]',I*J,1) .*
                                      reshape(Vn1[:,:,indx_k],I*J,N-1), dims=2)
        			qk 			= Delta*uk_stacked + Vk_stacked + Delta*Vkp_stacked
	        		Vn2k_stacked = Bik_all[kk]*qk
        			Vn2[:,:,kk] = reshape(Vn2k_stacked,I,J,1)
		        end
	        	VHIS_delta = Vn2 - Vn1
        		Vn1 		= Vn2
        		dist 		= maximum(abs.(VHIS_delta))

                if dist < crit_HIS
		        	break
        		end
	        end
        end


        # Check for convergence
        V_Delta = Vn1 - Vn
        Vn = Vn1
        dist = maximum(abs.(V_Delta))#maximum(maximum(maximum(abs.(V_Delta))))
        @show dist

        if dist < crit_HJB
 	        println("Value Function Converged, Iteration = ", nn)
	        break
        end
    end

	# Store value function
	V_0 = Vn

	#----------------------------------------------------------------
	# Compute stationary distribution and update aggregate KL
	#----------------------------------------------------------------

	#------------------------------------------------------------
	# Find new stationary distribution associated with decision rules
	#------------------------------------------------------------
	A   = aau + bbu
	λ0 	= lambda - diagm(0 => diag(lambda))  		# transition matrix with diagonal killed
	λ0p = λ0'

    gg_tilde = dab_g_tilde_mat * gg
    gg1      = Array{Float64}(undef, I_g*J_g,N)
    gg       = Array{Float64}(undef, I_g*J_g,N)
    g        = Array{Float64,3}(undef, I_g,J_g,N)

    K_supply  = 0.
    L_supply  = 0.
    KL_supply = 0.
    KL        = 0.

    # Iterate
    for nn = 1:maxit_KFE

        gg_tilde = dab_g_tilde_mat * gg
        gg1             = Array{Float64}(undef, I_g*J_g,N)

        for kk = 1:N
            Ak = A[1+(kk-1)*(I_g*J_g):kk*(I_g*J_g),1+(kk-1)*(I_g*J_g):kk*(I_g*J_g)]
            death_inflow = zeros(I_g, J_g)
            death_inflow[b.==0, 1] .= sum(gg_tilde[1+(kk-1)*(I_g*J_g):kk*(I_g*J_g)])
            death_inflow = reshape(death_inflow,I_g*J_g,1)
            gk_sum = sum(repeat(λ0p[kk,:]',I_g*J_g,1) .* reshape(gg_tilde,I_g*J_g,N), dims=2)
            gg1[:,kk] = (my_speye(I_g*J_g) - Delta_KFE * Ak' - Delta_KFE * (lambda[kk,kk] - ddeath) *
                         my_speye(I_g*J_g))\(gg_tilde[1+(kk-1)*(I_g*J_g):kk*(I_g*J_g)] +
                                             Delta_KFE*gk_sum + Delta_KFE*ddeath*death_inflow)
        end

           gg1             = reshape(gg1,I_g*J_g*N,1)
           gg1_sum         = sum(gg1)
           gg1             = gg1 ./ gg1_sum
           gg1             = dab_g_tilde_mat \ gg1

           dist = maximum(abs.(gg1-gg))
           if mod(nn,100) == 0
               println("Iteration ",nn,", distance is ", dist)
           end
           if dist < crit_KFE
                println("Distribution Converged, Iteration = ", nn)
               break
           end
           gg     = gg1

        g    = reshape(gg,I_g,J_g,N)
        end

        #------------------------------------------------------------
        # Update guess of KL
        #------------------------------------------------------------

        # Capital supply
        if K_liquid == 1
            K_supply = sum(g.*(a_g_grid + b_g_grid).*dab_tilde_grid)
        else
            K_supply = sum(g.*a_g_grid.*dab_g_tilde_grid)
        end

        # Labor supply
        L_supply = sum(y_g_grid .* l_g_grid .* g .* dab_g_tilde_grid)

        # Capital-Labor ratio
        KL_supply = K_supply / L_supply

        # Check for convergence and update
        gap = KL - KL_supply
        println("The current gap is ", gap)

        if abs(gap) > crit_KL
            KL = relax_KL*KL_supply + (1-relax_KL)*KL
        else
            println("I have found the steady state, Iteration = ", ii)
            break
        end
    end


    #compute_savings()
    A = aau + bbu + ccu
    dab_small = reshape(dab_g_tilde_grid[:,:,1], I_g*J_g, 1)
    loc = findall(!iszero, b .== 0)
    dab_small      = dab_small ./ dab_small[loc] * ddeath
    dab_small[loc] .= 0.0
    death_process        = -ddeath * my_speye(I_g*J_g)
    death_process[loc,:] = vec(dab_small)
    death_process        = kron(my_speye(N),death_process)

    g_new = A' * gg_tilde
    g_new = dab_g_tilde_mat \ g_new
    g_new = g_new + death_process * gg
    g_new = reshape(g_new ,I_g, J_g, N)
    g_new = sum(sum(g_new, dims=1), dims=3)
    save_a = dot(g_new, a_g) #g_new * a_g

    g_new = A' * gg_tilde
    g_new = dab_g_tilde_mat\g_new
    g_new = g_new + death_process * gg
    g_new = reshape(g_new, I_g, J_g, N)
    g_new = sum(sum(g_new, dims=2), dims=3)
    save_b = dot(vec(g_new), b_g) # g_new' * b_g
    #compute_savings()

    # Rename variables in steady state
    m[:KL_SS]			= KL
    m[:r_a_SS]			= aalpha * (m[:KL_SS] ^ (aalpha - 1)) - ddelta
    m[:w_SS]			= (1 - aalpha) * (m[:KL_SS] ^ aalpha)
    m[:K_SS]		    = m[:KL_SS] * L_supply
    m[:u_SS]			= u
    m[:c_SS]			= c
    m[:d_SS]			= d
    m[:V_SS]			= Vn
    m[:g_SS]			= g
    m[:dab]				= dab_tilde_grid
    m[:dab_g]           = dab_g_tilde_grid
    m[:C_SS]			= sum(g .* c_g .* dab_g_tilde_grid)
    m[:C_Var_SS]        = sum(g .* log.(c_g).^2 .* dab_g_tilde_grid) -
                                 sum(g .* log.(c_g) .* dab_g_tilde_grid)^2
    m[:I_SS]			= save_a
    m[:B_SS]			= sum(g .* b_g_grid .* dab_g_tilde_grid)
    m[:Y_SS]			= (K_supply ^ aalpha) * (L_supply ^ (1 - aalpha))
    m[:n_SS]			= L_supply
    m[:earn_SS]         = log.((1-tau_I) * w * l_g_grid .* y_g_grid + b_g_grid .*
                               (r_b_g_grid .+ ddeath*pam) .+ trans_g_grid .+ a_g_grid .*
                               (r_a_g_grid .+ ddeath*pam))
    m[:earn_Var_SS]     = sum(g .* m[:earn_SS].^2 .* dab_g_tilde_grid) -
                                sum(g .* m[:earn_SS] .* dab_g_tilde_grid)^2
    m[:IcF_SS]	 = IcF
    m[:IcB_SS]	 = IcB
    m[:Ic0_SS]	 = Ic0
    m[:IcFB_SS]	 = IcFB
    m[:IcBF_SS]	 = IcBF
    m[:IcBB_SS]	 = IcBB
    m[:Ic00_SS]	 = Ic00
    m[:a_g_0pos] = a_g_0pos[1]
    m[:b_g_0pos] = b_g_0pos[1]

    ###
    # Consumption by hand-to-mouth status
    ###
    WHTM_indicator      = zeros(I_g, J_g, N)
    WHTM_indicator[b_g_0pos:b_g_0pos+1, a_g_0pos+2:end,:] .= 1.0

    m[:WHTM_indicator] = WHTM_indicator
    m[:WHTM_SS]        = sum(g[:] .* WHTM_indicator[:] .* dab_g_tilde_grid[:])
    m[:C_WHTM_SS]      = sum(WHTM_indicator[:] .* c_g[:] .* g[:] .* dab_g_tilde_grid[:])

    PHTM_indicator      = zeros(I_g, J_g, N)
    PHTM_indicator[b_g_0pos:b_g_0pos+1,a_g_0pos:a_g_0pos+2:end,:] .= 1.0

    m[:PHTM_indicator] = PHTM_indicator
    m[:PHTM_SS]        = sum(g[:] .* PHTM_indicator[:] .* dab_g_tilde_grid[:])
    m[:C_PHTM_SS]      = sum(c_g[:] .* g[:] .* PHTM_indicator[:] .* dab_g_tilde_grid[:])

    NHTM_indicator      = zeros(I_g,J_g,N)
    NHTM_indicator[(b_g_0pos+3):end,2:end,:] .= 1.0

    m[:NHTM_indicator_SS] = NHTM_indicator
    m[:NHTM_SS]           = sum(g[:] .* NHTM_indicator[:] .* dab_g_tilde_grid[:])
    m[:C_NHTM_SS]         = sum(c_g[:] .* g[:] .* NHTM_indicator[:] .* dab_g_tilde_grid[:]) /
                        m[:NHTM_SS]
    m[:r_a_grid]	      = repeat([r_a], I_g, J_g, N)

    # Collect variables into vector
    vars_SS                      = zeros(n_v + n_g + n_p + n_Z, 1)
    vars_SS[1:n_v]               = reshape(m[:V_SS], I*J*N, 1)

    gg_SS                        = reshape(m[:g_SS], I_g*J_g*N, 1)
    m[:gg_SS] = gg_SS

    vars_SS[n_v+1:n_v+n_g,1]     = gg_SS[1:I_g*J_g*N-1]
    vars_SS[n_v+n_g+1,1]         = m[:K_SS]
    vars_SS[n_v+n_g+2,1]         = m[:r_b_SS]

    if aggregate_variables == 1
        vars_SS[n_v+n_g+3,1]         = m[:Y_SS]
        vars_SS[n_v+n_g+4,1]         = m[:C_SS]
    elseif distributional_variables == 1
        vars_SS[n_v+n_g+3,1]        = m[:C_Var_SS]
    	vars_SS[n_v+n_g+4,1]        = m[:earn_Var_SS]
    elseif distributional_variables_1 == 1
        vars_SS[n_v+n_g+3,1]        = m[:C_WHTM_SS]
    	vars_SS[n_v+n_g+4,1]        = m[:C_PHTM_SS]
        #vars_SS[n_v+n_g+3,1]        = C_NHTM_SS
    end
    vars_SS[n_v+n_g+n_p+1,1] = 0.0

    # Compute illiquid wedge
    m[:illiquid_wedge] = m[:r_a_SS] - m[:r_b_SS]

    # Save SS variables
    m[:var_SS] = var_SS

    return m
end
