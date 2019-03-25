"""
``
eqcond(m::TwoAssetHANK)
```

Expresses the equilibrium conditions in canonical form using Γ0, Γ1, C, Ψ, and Π matrices.

### Outputs

* `Γ0` (`n_states` x `n_states`) holds coefficients of current time states.
* `Γ1` (`n_states` x `n_states`) holds coefficients of lagged states.
* `C`  (`n_states` x `1`) is a vector of constants
* `Ψ`  (`n_states` x `n_shocks_exogenous`) holds coefficients of iid shocks.
* `Π`  (`n_states` x `n_states_expectational`) holds coefficients of expectational states.
"""
function eqcond(m::TwoAssetHANK)
    # Read in parameters
    aalpha = m[:aalpha].value
    ddelta = m[:ddelta].value
    ddeath = m[:ddeath].value
    rrho   = m[:rrho].value
    chi0   = m[:chi0].value
    chi1   = m[:chi1].value
    chi2   = m[:chi2].value
    a_lb   = m[:a_lb].value
    kappa  = m[:kappa].value
    pam    = m[:pam].value
    xxi    = m[:xxi].value
    ggamma = m[:ggamma].value
    tau_I  = m[:tau_I].value
    trans  = m[:trans].value
    dab    = m[:dab].value
    dab_g  = m[:dab_g].value
    n_SS   = m[:n_SS].value

    nnu_aggZ    = m[:nnu_aggZ].value
    ssigma_aggZ = m[:ssigma_aggZ].value

    # Set liquid rates
    r_b_SS       = m[:r_b_SS].value
    r_b_borr_SS  = m[:r_b_borr_SS].value
    borrwedge_SS = m[:borrwedge_SS].value

    lambda                     = get_setting(m, :lambda)
    K_liquid                   = get_setting(m, :K_liquid)
    aggregate_variables        = get_setting(m, :aggregate_variables)
    distributional_variables   = get_setting(m, :distributional_variables)
    distributional_variables_1 = get_setting(m, :distributional_variables_1)
    interp_decision            = get_setting(m, :interp_decision)
    permanent                  = get_setting(m, :permanent)

    I      = get_setting(m, :I)
    J      = get_setting(m, :J)
    a_g    = get_setting(m, :a_g)
    b_g    = get_setting(m, :b_g)
    I_g    = get_setting(m, :I_g)
    J_g    = get_setting(m, :J_g)
    N      = get_setting(m, :N)
    a      = get_setting(m, :a)
    b      = get_setting(m, :b)
    y      = get_setting(m, :y)
    y_dist = get_setting(m, :y_dist)
    y_mean = get_setting(m, :y_mean)
    KL     = get_setting(m, :KL_0)
    r_b_fix= get_setting(m, :r_b_fix)

    # Set liquid rates
    r_b      = r_b_SS
    r_b_borr = r_b_borr_SS

    # Compute prices associated with initial guess of KL
    w	= (1 - aalpha) * (KL ^ aalpha)
    r_a	= aalpha * (KL ^ (aalpha - 1)) - ddelta

a_grid, a_g_grid, b_grid, b_g_grid, y_grid, y_g_grid, r_a_grid, r_b_grid, r_a_g_grid, r_b_g_grid, daf_grid, daf_g_grid, dab_grid, dab_g_grid, dab_tilde_grid, dab_g_tilde_grid, dab_g_tilde_mat, dab_g_tilde, dbf_grid, dbf_g_grid, dbb_grid, dbb_g_grid, trans_grid, trans_g_grid, l_grid, l_g_grid, w_grid = set_grids(a, b, a_g, b_g, y, I, J, I_g, J_g, N, w, r_a, r_b, r_b_borr, trans)

    n_v = get_setting(m, :n_v)
    n_g = get_setting(m, :n_g)
    n_p = get_setting(m, :n_p)
    n_Z = get_setting(m, :n_Z)

    IcF_SS   = get_setting(m, :IcF_SS)
    IcB_SS   = get_setting(m, :IcB_SS)
    Ic0_SS   = get_setting(m, :Ic0_SS)
    IcFB_SS  = get_setting(m, :IcFB_SS)
    IcBF_SS  = get_setting(m, :IcBF_SS)
    IcBB_SS  = get_setting(m, :IcBB_SS)
    Ic00_SS  = get_setting(m, :Ic00_SS)
    nVars    = get_setting(m, :nVars)
    nEErrors = get_setting(m, :nEErrors)

    vars_SS = m[:vars_SS].value

    # Computation taken from inside get_residuals
    cc  = kron(lambda, my_speye(I*J))
    ccu = kron(lambda, my_speye(I_g*J_g))

    dab_aux   = reshape(dab,I*J*N,1)
    dab_g_aux = reshape(dab_g,I_g*J_g*N,1)

    # Which direction to use
    IcF = IcF_SS
    IcB = IcB_SS
    Ic0 = Ic0_SS

    # Which direction to use
    IcFB = IcFB_SS
    IcBF = IcBF_SS
    IcBB = IcBB_SS
    Ic00 = Ic00_SS

    @inline function get_residuals(vars::Vector{T}) where {T<:Real}
        # Unpack variables
        V         = reshape(vars[1:n_v] .+ vars_SS[1:n_v], I, J, N)  # value function
        g         = vars[n_v+1:n_v+n_g] .+ vars_SS[n_v+1:n_v+n_g]    # distribution
        g_end     = (1 - sum(g .* dab_g_aux[1:end-1])) / dab_g[I_g,J_g,N]
        gg        = [g;g_end]
        K         = vars[n_v+n_g+1] + vars_SS[n_v+n_g+1]    # aggregate capital
        r_b       = vars[n_v+n_g+2] + vars_SS[n_v+n_g+2]

        if aggregate_variables == 1
            aggY     = vars[n_v+n_g+3] + vars_SS[n_v+n_g+3] # aggregate output
            aggC     = vars[n_v+n_g+4] + vars_SS[n_v+n_g+4] # aggregate consumption
        elseif distributional_variables == 1
            C_Var    = vars[n_v+n_g+3] + vars_SS[n_v+n_g+3] # consumption inequality
            earn_Var = vars[n_v+n_g+4] + vars_SS[n_v+n_g+4] # earnings inequality
        elseif distributional_variables_1 == 1
            C_WHTM  = vars[n_v+n_g+3] + vars_SS[n_v+n_g+3] # consumption of wealthy hand-to-mouth
            C_PHTM  = vars[n_v+n_g+4] + vars_SS[n_v+n_g+4] # consumption of poor hand-to-mouth
        end

        aggZ       = vars[n_v+n_g+n_p+1] + vars_SS[n_v+n_g+n_p+1] # aggregate Z

        V_Dot      = vars[nVars+1:nVars+n_v]
        g_Dot      = vars[nVars+n_v+1:nVars+n_v+n_g]
        aggZ_Dot   = vars[nVars+n_v+n_g+n_p+1]

        VEErrors   = vars[2*nVars+1:2*nVars+n_v]
        aggZ_Shock = vars[2*nVars+nEErrors+1]

        # Prices
        if permanent == 0
            w   = (1 - aalpha) * (K ^ aalpha) * exp(aggZ) ^ (1-aalpha) * n_SS ^ (-aalpha)
            r_a = aalpha * (K ^ (aalpha - 1)) * (exp(aggZ) * n_SS) ^ (1 - aalpha) - ddelta
        elseif permanent == 1
            w   = (1 - aalpha) * (K ^ aalpha) * (n_SS) ^ (-aalpha)
            r_a = aalpha * (K ^ (aalpha - 1)) * (n_SS) ^ (1 - aalpha) - ddelta
        end

        # Auxiliary variables
        r_b_borr = r_b .+ borrwedge_SS

        # Other necessary objects
        r_a_grid     = repeat([r_a],I,J,N)

        l_grid       = permutedims(repeat(ones(N,1),1,I,J), [2 3 1])
        y_shock      = y .* exp.(kappa * aggZ * (y .- y_mean) ./ std(y))
        y_shock_mean = dot(y_shock, y_dist) #y_shock * y_dist
        y_shock      = y_shock ./ y_shock_mean .* y_mean
        y_grid       = permutedims(repeat(y_shock',1,I,J), [2 3 1])

        #----------------------------------------------------------------
        # HJB Equation
        #----------------------------------------------------------------
        # Preparations
        VaF   = similar(V) #0. * V
        VaB   = similar(V) #0. * V
        Vamin = 0.0

        # Forward difference
        VaF[:,1:J-1,:] = (V[:,2:J,:] - V[:,1:J-1,:]) ./ daf_grid[:,1:J-1,:]
        VaF[:,1:J-1,:] = max.(VaF[:,1:J-1,:], Vamin)
        VaF[:,J,:]    .= 0.0

        # Backward difference
        VaB[:,2:J,:] = (V[:,2:J,:] - V[:,1:J-1,:]) ./ dab_grid[:,2:J,:]
        VaB[:,2:J,:] = max.(VaB[:,2:J,:], Vamin)
        VaB[:,1,:]  .= 0.0

        # Preparations (necessary to ensure that everything is a dual number,
        # required by derivative software)
        VbF   = similar(V) #0 * V
        VbB   = similar(V) #0 * V
        Vbmin = 1e-8

        # Forward difference
        VbF[1:I-1,:,:] = (V[2:I,:,:] - V[1:I-1,:,:]) ./ dbf_grid[1:I-1,:,:]
        VbF[1:I-1,:,:] = max.(VbF[1:I-1,:,:], Vbmin)
        VbF[I,:,:]    .= 0.0

        # Backward difference
        VbB[2:I,:,:] = (V[2:I,:,:] - V[1:I-1,:,:]) ./ dbb_grid[2:I,:,:]
        VbB[2:I,:,:] = max.(VbB[2:I,:,:], Vbmin)
        VbB[1,:,:]  .= 0.0

        #----------------------------------------------------------------
        # Consumption decision
        #----------------------------------------------------------------

        # Preparations
        cF  = similar(V) #0 * V
        sF  = similar(V) #0 * V
        HcF = similar(V) #0 * V

        cB  = similar(V) #0 * V
        sB  = similar(V) #0 * V
        HcB = similar(V) #0 * V

        c0  = similar(V) #0 * V
        #s0  = similar(V) #0 * V
        Hc0 = similar(V) #0 * V

        # Decisions conditional on a particular direction of derivative
        cF[1:I-1,:,:] = VbF[1:I-1,:,:] .^ (-1/ggamma)
        cF[I,:,:]    .= 0.0 #zeros(1,J,N)

        if permanent == 1
            sF[1:I-1,:,:] = real(((1-xxi)-tau_I) * w * l_grid[1:I-1,:,:] .* y_grid[1:I-1,:,:] .+
                                 b_grid[1:I-1,:,:].*(r_b_grid[1:I-1,:,:] .+ ddeath*pam - aggZ) .+
                                 trans_grid[1:I-1,:,:] - cF[1:I-1,:,:])
        elseif permanent == 0
            sF[1:I-1,:,:] = real(((1-xxi)-tau_I) * w * l_grid[1:I-1,:,:] .* y_grid[1:I-1,:,:] .+
                                 b_grid[1:I-1,:,:] .* (r_b_grid[1:I-1,:,:] .+ ddeath*pam) .+
                                 trans_grid[1:I-1,:,:] - cF[1:I-1,:,:])
        end

        sF[I,:,:]     .= 0.0 #zeros(1,J,N)
        HcF[1:I-1,:,:] = u_fn(cF[1:I-1,:,:], ggamma) .+ VbF[1:I-1,:,:] .* sF[1:I-1,:,:]
        HcF[I,:,:]    .= -1e12 #*ones(1,J,N)
        validF         = (sF .> 0)

        cB[2:I,:,:]    = VbB[2:I,:,:].^(-1/ggamma)

        if permanent == 1
            cB[1,:,:] = ((1-xxi)-tau_I) * w * l_grid[1,:,:] .* y_grid[1,:,:] .+
                b_grid[1,:,:] .* (r_b_grid[1,:,:] .+ ddeath*pam - aggZ) .+ trans_grid[1,:,:]
            sB[2:I,:,:]     = ((1-xxi)-tau_I) * w * l_grid[2:I,:,:] .* y_grid[2:I,:,:] .+
                b_grid[2:I,:,:] .* (r_b_grid[2:I,:,:] .+ ddeath*pam - aggZ) .+
                trans_grid[2:I,:,:] - cB[2:I,:,:]
        elseif permanent == 0
            cB[1,:,:]     = real(((1-xxi)-tau_I) * w * l_grid[1,:,:] .* y_grid[1,:,:] .+
                b_grid[1,:,:] .* (r_b_grid[1,:,:] .+ ddeath*pam) .+ trans_grid[1,:,:])
            sB[2:I,:,:] = real(((1-xxi)-tau_I) * w * l_grid[2:I,:,:] .*
                y_grid[2:I,:,:] .+ b_grid[2:I,:,:] .* (r_b_grid[2:I,:,:] .+ ddeath*pam) .+
                trans_grid[2:I,:,:] - cB[2:I,:,:])
        end

        sB[1,:,:] .= 0.0 #zeros(1,J,N)
        HcB        = u_fn(cB, ggamma) .+ VbB .* sB
        validB     = (sB .< 0)

        if permanent == 1
            c0 = real(((1-xxi)-tau_I) * w * l_grid .* y_grid .+
                b_grid .* (r_b_grid .+ ddeath*pam - aggZ) .+ trans_grid)
        elseif permanent == 0
            c0 = real(((1-xxi)-tau_I) * w * l_grid .* y_grid .+
                b_grid .* (r_b_grid .+ ddeath*pam) .+ trans_grid)
        end

        #s0  = zeros(I,J,N)
        Hc0 = u_fn(c0, ggamma)


        # Optimal consumption and liquid savings
        c = IcF .* cF .+ IcB .* cB .+ Ic0 .* c0
        s = IcF .* sF .+ IcB .* sB #.+ Ic0 .* s0
        u = u_fn(c, ggamma)

        #----------------------------------------------------------------
        # Deposit decision
        #----------------------------------------------------------------

        # Preparations
        dFB  = similar(V) #0*V
        HdFB = similar(V) #0*V
        dBF  = similar(V) #0*V
        HdBF = similar(V) #0*V
        dBB  = similar(V) #0*V
        HdBB = similar(V) #0*V

        # Decisions conditional on a particular direction of derivative
        dFB[2:I,1:J-1,:]  = opt_deposits(VaF[2:I,1:J-1,:], VbB[2:I,1:J-1,:],
                                         a_grid[2:I,1:J-1,:], chi0, chi1, chi2, a_lb)
        dFB[:,J,:]        .= 0.0 #zeros(I,1,N)
        dFB[1,1:J-1,:]    .= 0.0 #zeros(1,J-1,N)
        HdFB[2:I,1:J-1,:]  = VaF[2:I,1:J-1,:] .* dFB[2:I,1:J-1,:] - VbB[2:I,1:J-1,:] .*
            (dFB[2:I,1:J-1,:] .+ adj_cost_fn(dFB[2:I,1:J-1,:],a_grid[2:I,1:J-1,:],
                                             chi0, chi1, chi2, a_lb))
        HdFB[:,J,:]     .= -1.0e12 #* ones(I,1,N)
        HdFB[1,1:J-1,:] .= -1.0e12 #* ones(1,J-1,N)
        validFB          = (dFB .> 0) .* (HdFB .> 0)

        dBF[1:I-1,2:J,:]  = opt_deposits(VaB[1:I-1,2:J,:], VbF[1:I-1,2:J,:],
                                            a_grid[1:I-1,2:J,:], chi0, chi1, chi2, a_lb)
        dBF[:,1,:]        .= 0.0 #zeros(I,1,N)
        dBF[I,2:J,:]      .= 0.0 #zeros(1,J-1,N)
        HdBF[1:I-1,2:J,:]  = VaB[1:I-1,2:J,:] .* dBF[1:I-1,2:J,:] - VbF[1:I-1,2:J,:] .*
            (dBF[1:I-1,2:J,:] .+ adj_cost_fn(dBF[1:I-1,2:J,:], a_grid[1:I-1,2:J,:],
                                             chi0, chi1, chi2, a_lb))
        HdBF[:,1,:]   .= -1.0e12 #* ones(I,1,N)
        HdBF[I,2:J,:] .= -1.0e12 #* ones(1,J-1,N)
        validBF        = (dBF .<= -adj_cost_fn(dBF, a_grid, chi0, chi1, chi2, a_lb)) .*
            (HdBF .> 0)

        VbB[1,2:J,:]  = u_fn(cB[1,2:J,:], ggamma)
        dBB[:,2:J,:]  = opt_deposits(VaB[:,2:J,:], VbB[:,2:J,:], a_grid[:,2:J,:],
                                           chi0, chi1, chi2, a_lb)
        dBB[:,1,:]   .= 0.0 #zeros(I,1,N)
        HdBB[:,2:J,:] = VaB[:,2:J,:] .* dBB[:,2:J,:] - VbB[:,2:J,:] .*
            (dBB[:,2:J,:] .+ adj_cost_fn(dBB[:,2:J,:], a_grid[:,2:J,:], chi0, chi1, chi2, a_lb))
        HdBB[:,1,:]   .= -1.0e12 #* ones(I,1,N)
        validBB       = (dBB .> -adj_cost_fn(dBB, a_grid, chi0, chi1, chi2, a_lb)) .*
            (dBB .<= 0) .* (HdBB .> 0)

        # Optimal deposit decision
        d = IcFB .* dFB .+ IcBF .* dBF .+ IcBB .* dBB #.+ Ic00 .* zeros(I,J,N)

        ## Interpolate
        d_g = reshape(interp_decision * vec(d), I_g, J_g, N)
        s_g = reshape(interp_decision * vec(s), I_g, J_g, N)
        c_g = reshape(interp_decision * vec(c), I_g, J_g, N)

        # Compute drifts for KFE
        @time audriftB, budriftB, audriftF, budriftF, adriftB, bdriftB, adriftF, bdriftF = catch_my_drifts(I_g, permanent, ddeath, pam, xxi, d_g, a_g_grid, r_a_g_grid, w,
                                 l_g_grid, y_g_grid, s_g, chi0, chi1, chi2,
                                 a_lb, a_grid, r_a_grid, l_grid, y_grid, aggZ, d, s)

        # Derive transition matrices
        @time aa, bb, aau, bbu = transition_deriva(I_g, J_g, N, I, J, ddeath, pam, xxi, w,
                                                   chi0, chi1,
                                             chi2, a_lb, l_grid, l_g_grid, y_grid, y_g_grid, d,
                                             dab_grid, daf_grid, dab_g_grid, daf_g_grid, dbb_grid,
                                             dbf_grid, dbb_g_grid, dbf_g_grid, d_g, a_grid,
                                             a_g_grid, s, s_g, r_a_grid, r_a_g_grid,
                                             audriftB, budriftB, audriftF, budriftF, adriftB,
                                             bdriftB, adriftF, bdriftF)
        # full transition matrix
        A = aa + bb + cc

        #----------------------------------------------------------------
        # KFE
        #----------------------------------------------------------------

        # KFE
        AT = aau + bbu + ccu
        AT = AT'
        gg_tilde            = dab_g_tilde_mat * gg
        gIntermediate       = AT * gg_tilde
        dab_g_tilde_mat_inv = spdiagm(0 => vec(repeat(1.0 ./ dab_g_tilde, N, 1)))
        gIntermediate       = dab_g_tilde_mat_inv * gIntermediate

        dab_g_small = reshape(dab_g[:,:,1],I_g*J_g,1)

        loc = findall(b .== 0)
        dab_g_small           = dab_g_small ./ dab_g_small[loc] * ddeath
        dab_g_small[loc]     .= 0.0
        death_process         = -ddeath * my_speye(I_g * J_g)
        death_process[loc,:]  = vec(dab_g_small)
        death_process         = kron(my_speye(N), death_process)

        gIntermediate = gIntermediate + death_process * gg

        # find death-corrected savings
        a_g_grid_aux = reshape(a_g_grid,I_g*J_g*N,1)
        a_save       = sum(gIntermediate .* dab_g_aux .* a_g_grid_aux)

        b_g_grid_aux = reshape(b_g_grid,I_g*J_g*N,1)
        b_save       = sum(gIntermediate .* dab_g_aux .* b_g_grid_aux)

        # consumption for low types and high types
        c_low     = c[:,:,1]
        dab_low   = dab[:,:,1]
        c_high    = c[:,:,2]
        dab_high  = dab[:,:,2]

        #----------------------------------------------------------------
        # Compute equilibrium conditions
        #----------------------------------------------------------------
        # HJB equation
        if permanent == 0
            hjbResidual = reshape(u, I*J*N, 1) + A * reshape(V, I*J*N, 1) + V_Dot + VEErrors -
                (rrho + ddeath) * reshape(V,I*J*N,1)
        elseif permanent == 1
            hjbResidual = reshape(u, I*J*N, 1) + A * reshape(V, I*J*N, 1) + V_Dot + VEErrors -
                (rrho + ddeath - (1 - ggamma) * aggZ) * reshape(V,I*J*N,1)
        end

        # KFE
        gResidual = g_Dot - gIntermediate[1:n_g, 1]

        K_out = 0.0
        if K_liquid == 1
            K_out = sum((vec(a_g_grid) .+ vec(b_g_grid)) .* gg .* vec(dab_g))
        else
            K_out = sum(vec(a_g_grid) .* gg .* vec(dab_g))
        end

        r_b_out = 0.0
        if r_b_fix      == 1
            r_b_out          = r_b_SS
        elseif r_b_phi  == 1
            r_b_out          = sum(vec(b_g_grid) .* gg .* vec(dab_g))
        elseif B_fix    == 1
            r_b_out          = 0.0
        elseif K_liquid == 1
            r_b_out          = r_a_out - illiquid_wedge
        end

        if aggregate_variables == 1
                aggY_out = (K ^ aalpha) * (n_SS ^ (1 - aalpha))
                aggC_out = sum(vec(c_g) .* gg .* vec(dab_g))
        elseif distributional_variables == 1
                C_Var_out = sum(log(vec(c_g)).^2 .* gg .* vec(dab_g)) -
                    sum(log(vec(c_g)) .* gg .* vec(dab_g)) ^ 2
                earn = log((1-tau_I) * w * l_g_grid .* y_g_grid + b_g_grid .*
                           (r_b_g_grid + ddeath*pam) + trans_grid + a_g_grid .*
                           (r_a_g_grid + ddeath*pam))
                earn_Var_out = sum(vec(earn).^2 .* gg .* vec(dab_g)) -
                    sum(vec(earn) .* gg .* vec(dab_g)) ^ 2
        elseif distributional_variables_1 == 1
                WHTM_indicator      = zeros(I_g,J_g,N)
                WHTM_indicator[b_g_0pos:b_g_0pos+1,a_g_0pos+2:end,:] .= 1.
                WHTM_out            = sum(vec(WHTM_indicator) .* gg .* vec(dab_g))
                C_WHTM_out          = sum(vec(WHTM_indicator) .* vec(c_g) .* gg .* vec(dab_g))

            if permanent == 0
                PHTM_indicator      = zeros(I_g,J_g,N)
                PHTM_indicator[b_g_0pos:b_g_0pos+1,a_g_0pos:a_g_0pos+2:end,:] .= 1.
                PHTM_out            = sum(vec(PHTM_indicator) .* gg .* vec(dab_g))
                C_PHTM_out          = sum(vec(PHTM_indicator) .* vec(c_g) .* gg .* vec(dab_g))

                NHTM_indicator      = zeros(I_g,J_g,N)
                NHTM_indicator[b_g_0pos+2:end,a_g_0pos+2:end,:] .= 1.
                NHTM_indicator      = zeros(I_g,J_g,N)
                NHTM_indicator[b_g_0pos+3:end,2:end,:] .= 1.
                NHTM_out            = sum(vec(NHTM_indicator) .* gg .* vec(dab_g))
                C_NHTM_out          = sum(vec(NHTM_indicator) .* vec(c_g) .* gg .* vec(dab_g))/NHTM_out

            elseif permanent == 1
                PHTM_indicator      = zeros(I_g,J_g,N)
                PHTM_indicator[b_g_0pos:b_g_0pos+1,a_g_0pos:a_g_0pos+2:end,:] .= 1.
                PHTM_out            = sum(vec(PHTM_indicator) .* gg .* vec(dab_g))
                C_PHTM_out          = sum(vec(PHTM_indicator) .* vec(c_g) .* gg .* vec(dab_g))

                NHTM_indicator      = zeros(I_g,J_g,N)
                NHTM_indicator[b_g_0pos+2:end,:,:] .= 1.
                NHTM_out            = sum(vec(NHTM_indicator) .* gg .* vec(dab_g))
                C_NHTM_out          = sum(vec(NHTM_indicator) .* vec(c_g) .* gg .* vec(dab_g))/NHTM_out
            end
        end

        K_Residual = K_out - K
        if r_b_fix == 1
            r_b_Residual = r_b_out - r_b
        elseif r_b_phi == 1
            r_b_Residual = r_b_out - B_SS * exp(1/pphi * (r_b - r_b_SS))
        elseif B_fix == 1
            r_b_Residual = r_b_out - b_save
        elseif K_liquid == 1
            r_b_Residual = r_b_out - r_b
        end

        Y_Residual        = Array{Float64}(undef, 0)
        C_Residual        = Array{Float64}(undef, 0)

        C_Var_Residual    = Array{Float64}(undef, 0)
        earn_Var_Residual = Array{Float64}(undef, 0)

        C_WHTM_Residual   = Array{Float64}(undef, 0)
        C_PHTM_Residual   = Array{Float64}(undef, 0)

        if aggregate_variables == 1
            Y_Residual          = aggY_out - aggY
            C_Residual          = aggC_out - aggC
        elseif distributional_variables == 1
            C_Var_Residual      = C_Var_out - C_Var
            earn_Var_Residual   = earn_Var_out - earn_Var
        elseif distributional_variables_1 == 1
            C_WHTM_Residual     = C_WHTM_out - C_WHTM
            C_PHTM_Residual     = C_PHTM_out - C_PHTM
        end

        # Law of motion for aggregate tfp shock
        aggZ_Residual = aggZ_Dot - (-nnu_aggZ * aggZ + ssigma_aggZ * aggZ_Shock)
        vResidual     = Array{Float64}(undef, 0)

        # Return equilibrium conditions
        if aggregate_variables == 1

            vResidual = [hjbResidual; gResidual; K_Residual; r_b_Residual; Y_Residual;
                         C_Residual; aggZ_Residual]
        elseif distributional_variables == 1

            vResidual = [hjbResidual; gResidual; K_Residual; r_b_Residual; C_Var_Residual;
                         earn_Var_Residual; aggZ_Residual]
        elseif distributional_variables_1 == 1

            vResidual = [hjbResidual; gResidual; K_Residual; r_b_Residual; C_WHTM_Residual;
                         C_PHTM_Residual; aggZ_Residual]
        else

            vResidual = [hjbResidual; gResidual; K_Residual; r_b_Residual; aggZ_Residual]
        end
        return vResidual
    end
#=
    function get_residuals



    # equilibrium conditions
    if aggregate_variables == 1
    vResidual = [hjbResidual; gResidual; K_Residual; r_b_Residual; Y_Residual;
                 C_Residual; aggZ_Residual]
    elseif distributional_variables == 1
        vResidual = [hjbResidual; gResidual; K_Residual; r_b_Residual; C_Var_Residual;
                     earn_Var_Residual; aggZ_Residual]
    elseif distributional_variables_1 == 1
        vResidual = [hjbResidual; gResidual; K_Residual; r_b_Residual; C_WHTM_Residual;
                     C_PHTM_Residual; aggZ_Residual]
    else
        vResidual = [hjbResidual; gResidual; K_Residual; r_b_Residual; aggZ_Residual]
    end
=#
    #####
    #nstates = n_states(m)
    #n_s_exp = n_shocks_expectational(m)
    #n_s_exo = n_shocks_exogenous(m)
    #vars = zeros(Float64, 2 * nstates + n_s_exp + n_s_exo)
    ####
    @time get_residuals(zeros(Float64, 2 * nVars + nEErrors + 1))
error()
    derivs = ForwardDiff.jacobian(get_residuals, zeros(Float64, 2 * nVars + nEErrors + 1))

    Γ1 = -derivs[:, 1:nstates]
    Γ0 =  derivs[:,   nstates +           1:2 * nstates]
    Π  = -derivs[:, 2*nstates +           1:2 * nstates + n_s_exp]
    Ψ  = -derivs[:, 2*nstates + n_s_exp + 1:2 * nstates + n_s_exp + n_s_exo]
    C  = zeros(nstates)

    if typeof(Ψ) == Vector{Float64}
        Ψ = reshape(Ψ, length(Ψ), n_s_exo)
    end

    return Γ0, Γ1, Ψ, Π, C
end
