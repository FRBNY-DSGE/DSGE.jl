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
    w2	= (1 - aalpha) * (KL ^ aalpha)
    r_a2	= aalpha * (KL ^ (aalpha - 1)) - ddelta

a_grid, a_g_grid, b_grid, b_g_grid, y_grid, y_g_grid, r_a_grid, r_b_grid, r_a_g_grid, r_b_g_grid, daf_grid, daf_g_grid, dab_grid, dab_g_grid, dab_tilde_grid, dab_g_tilde_grid, dab_g_tilde_mat, dab_g_tilde, dbf_grid, dbf_g_grid, dbb_grid, dbb_g_grid, trans_grid, trans_g_grid, l_grid, l_g_grid, w_grid = set_grids(a, b, a_g, b_g, y, I, J, I_g, J_g, N, w2, r_a2, r_b, r_b_borr, trans)

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
#=
    @inline function get_residuals_agg(vars::Vector{T}) where {T<:Real}
        # Unpack variables
        V         = reshape(vars[1:n_v] .+ vars_SS[1:n_v], I, J, N)  # value function
        g         = vars[n_v+1:n_v+n_g] .+ vars_SS[n_v+1:n_v+n_g]    # distribution
        g_end     = (1 - sum(g .* dab_g_aux[1:end-1])) / dab_g[I_g,J_g,N]
        gg        = [g;g_end]
        K         = vars[n_v+n_g+1] + vars_SS[n_v+n_g+1]    # aggregate capital
        r_b       = vars[n_v+n_g+2] + vars_SS[n_v+n_g+2]

        aggY     = vars[n_v+n_g+3] + vars_SS[n_v+n_g+3] # aggregate output
        aggC     = vars[n_v+n_g+4] + vars_SS[n_v+n_g+4] # aggregate consumption
        aggZ     = vars[n_v+n_g+n_p+1] + vars_SS[n_v+n_g+n_p+1] # aggregate Z

        V_Dot      = vars[nVars+1:nVars+n_v]
        g_Dot      = vars[nVars+n_v+1:nVars+n_v+n_g]
        aggZ_Dot   = vars[nVars+n_v+n_g+n_p+1]

        VEErrors   = vars[2*nVars+1:2*nVars+n_v]
        aggZ_Shock = vars[2*nVars+nEErrors+1]

        # Prices
        w   = (1 - aalpha) * (K ^ aalpha) * n_SS ^ (-aalpha) *
            ((permanent == 0) ? exp(aggZ) ^ (1-aalpha) : 1.)
        r_a = aalpha * (K ^ (aalpha - 1)) * (((permanent == 0) ? exp(aggZ) : 1.) *
                                             n_SS) ^ (1 - aalpha) - ddelta
        # Auxiliary variables
        r_b_borr = r_b .+ borrwedge_SS

        # Other necessary objects
        r_a_grid     = repeat([r_a],I,J,N)

        l_grid       = permutedims(repeat(ones(N,1),1,I,J), [2 3 1])
        y_shock      = y .* exp.(kappa * aggZ * (y .- y_mean) ./ std(y))
        y_shock_mean = dot(y_shock, y_dist) #y_shock * y_dist
        y_shock      = y_shock ./ y_shock_mean .* y_mean
        y_grid       = permutedims(repeat(y_shock',1,I,J), [2 3 1])

        # ripped out
        @time c, s, u, d, d_g, s_g, c_g = eqcond_helper(V, I, J, I_g, J_g, N, chi0, chi1, chi2,
                                                        a_lb, ggamma, permanent, interp_decision,
                                                        ddeath, pam, aggZ, xxi, tau_I, w, l_grid,
                                                        y_grid, b_grid, r_b_grid, trans_grid,
                                                        a_grid,
                                                        daf_grid, dab_grid, dbf_grid, dbb_grid,
                                                        IcF, IcB, Ic0, IcFB, IcBF, IcBB, Ic00)

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
        perm_mult   = (permanent == 0) ? rrho + ddeath : rrho + ddeath - (1 - ggamma) * aggZ
        hjbResidual = reshape(u, I*J*N, 1) + A * reshape(V, I*J*N, 1) + V_Dot + VEErrors -
            perm_mult * reshape(V,I*J*N,1)

        # KFE
        gResidual = g_Dot - gIntermediate[1:n_g, 1]

        K_out = 0.0
        if K_liquid == 1
            K_out = sum((vec(a_g_grid) .+ vec(b_g_grid)) .* gg .* vec(dab_g))
        else
            K_out = sum(vec(a_g_grid) .* gg .* vec(dab_g))
        end

        r_b_out = 0.0
        r_b_Residual = 0.0
        if r_b_fix      == 1
            r_b_out      = r_b_SS
            r_b_Residual = r_b_out - r_b
        elseif r_b_phi  == 1
            r_b_out      = sum(vec(b_g_grid) .* gg .* vec(dab_g))
            r_b_Residual = r_b_out - B_SS * exp(1/pphi * (r_b - r_b_SS))
        elseif B_fix    == 1
            r_b_out      = 0.0
            r_b_Residual = r_b_out - b_save
        elseif K_liquid == 1
            r_b_out      = r_a_out - illiquid_wedge
            r_b_Residual = r_b_out - r_b
        end

        aggY_out = (K ^ aalpha) * (n_SS ^ (1 - aalpha))
        aggC_out = sum(vec(c_g) .* gg .* vec(dab_g))
        K_Residual = K_out - K

        Y_Residual          = aggY_out - aggY
        C_Residual          = aggC_out - aggC

        # Law of motion for aggregate tfp shock
        aggZ_Residual = aggZ_Dot - (-nnu_aggZ * aggZ + ssigma_aggZ * aggZ_Shock)

        # Return equilibrium conditions
        return [hjbResidual; gResidual; K_Residual; r_b_Residual; Y_Residual; C_Residual;
                aggZ_Residual]
    end

    @inline function get_residuals_dist(vars::Vector{T}) where {T<:Real}
        # Unpack variables
        V         = reshape(vars[1:n_v] .+ vars_SS[1:n_v], I, J, N)  # value function
        g         = vars[n_v+1:n_v+n_g] .+ vars_SS[n_v+1:n_v+n_g]    # distribution
        g_end     = (1 - sum(g .* dab_g_aux[1:end-1])) / dab_g[I_g,J_g,N]
        gg        = [g;g_end]
        K         = vars[n_v+n_g+1] + vars_SS[n_v+n_g+1]    # aggregate capital
        r_b       = vars[n_v+n_g+2] + vars_SS[n_v+n_g+2]

        C_Var    = vars[n_v+n_g+3] + vars_SS[n_v+n_g+3] # consumption inequality
        earn_Var = vars[n_v+n_g+4] + vars_SS[n_v+n_g+4] # earnings inequality

        aggZ       = vars[n_v+n_g+n_p+1] + vars_SS[n_v+n_g+n_p+1] # aggregate Z

        V_Dot      = vars[nVars+1:nVars+n_v]
        g_Dot      = vars[nVars+n_v+1:nVars+n_v+n_g]
        aggZ_Dot   = vars[nVars+n_v+n_g+n_p+1]

        VEErrors   = vars[2*nVars+1:2*nVars+n_v]
        aggZ_Shock = vars[2*nVars+nEErrors+1]

        # Prices
        w   = (1 - aalpha) * (K ^ aalpha) * n_SS ^ (-aalpha) *
            ((permanent == 0) ? exp(aggZ) ^ (1-aalpha) : 1.)
        r_a = aalpha * (K ^ (aalpha - 1)) * (((permanent == 0) ? exp(aggZ) : 1.) *
                                             n_SS) ^ (1 - aalpha) - ddelta
        # Auxiliary variables
        r_b_borr = r_b .+ borrwedge_SS

        # Other necessary objects
        r_a_grid     = repeat([r_a],I,J,N)

        l_grid       = permutedims(repeat(ones(N,1),1,I,J), [2 3 1])
        y_shock      = y .* exp.(kappa * aggZ * (y .- y_mean) ./ std(y))
        y_shock_mean = dot(y_shock, y_dist) #y_shock * y_dist
        y_shock      = y_shock ./ y_shock_mean .* y_mean
        y_grid       = permutedims(repeat(y_shock',1,I,J), [2 3 1])

        # ripped out
        @time c, s, u, d, d_g, s_g, c_g = eqcond_helper(V, I, J, I_g, J_g, N, chi0, chi1, chi2,
                                                        a_lb, ggamma, permanent, interp_decision,
                                                        ddeath, pam, aggZ, xxi, tau_I, w, l_grid,
                                                        y_grid, b_grid, r_b_grid, trans_grid,
                                                        a_grid,
                                                        daf_grid, dab_grid, dbf_grid, dbb_grid,
                                                        IcF, IcB, Ic0, IcFB, IcBF, IcBB, Ic00)

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
        perm_mult = (permanent == 0) ? rrho + ddeath : rrho + ddeath - (1 - ggamma) * aggZ
        hjbResidual = reshape(u, I*J*N, 1) + A * reshape(V, I*J*N, 1) + V_Dot + VEErrors -
            perm_mult * reshape(V,I*J*N,1)

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
            r_b_Residual = r_b_out - r_b
        elseif r_b_phi  == 1
            r_b_out      = sum(vec(b_g_grid) .* gg .* vec(dab_g))
            r_b_Residual = r_b_out - B_SS * exp(1/pphi * (r_b - r_b_SS))
        elseif B_fix    == 1
            r_b_out          = 0.0
            r_b_Residual = r_b_out - b_save
        elseif K_liquid == 1
            r_b_out          = r_a_out - illiquid_wedge
            r_b_Residual = r_b_out - r_b
        end

        C_Var_out = sum(log(vec(c_g)).^2 .* gg .* vec(dab_g)) -
            sum(log(vec(c_g)) .* gg .* vec(dab_g)) ^ 2
        earn = log((1-tau_I) * w * l_g_grid .* y_g_grid + b_g_grid .*
                   (r_b_g_grid + ddeath*pam) + trans_grid + a_g_grid .*
                   (r_a_g_grid + ddeath*pam))
        earn_Var_out = sum(vec(earn).^2 .* gg .* vec(dab_g)) -
            sum(vec(earn) .* gg .* vec(dab_g)) ^ 2

        K_Residual = K_out - K

        C_Var_Residual      = C_Var_out - C_Var
        earn_Var_Residual   = earn_Var_out - earn_Var

        # Law of motion for aggregate tfp shock
        aggZ_Residual = aggZ_Dot - (-nnu_aggZ * aggZ + ssigma_aggZ * aggZ_Shock)

        # Return equilibrium conditions
        return [hjbResidual; gResidual; K_Residual; r_b_Residual; C_Var_Residual;
                earn_Var_Residual; aggZ_Residual]
    end

    @inline function get_residuals_dist1(vars::Vector{T}) where {T<:Real}
        # Unpack variables
        V         = reshape(vars[1:n_v] .+ vars_SS[1:n_v], I, J, N)  # value function
        g         = vars[n_v+1:n_v+n_g] .+ vars_SS[n_v+1:n_v+n_g]    # distribution
        g_end     = (1 - sum(g .* dab_g_aux[1:end-1])) / dab_g[I_g,J_g,N]
        gg        = [g;g_end]
        K         = vars[n_v+n_g+1] + vars_SS[n_v+n_g+1]    # aggregate capital
        r_b       = vars[n_v+n_g+2] + vars_SS[n_v+n_g+2]
        C_WHTM    = vars[n_v+n_g+3] + vars_SS[n_v+n_g+3] # consumption of wealthy hand-to-mouth
        C_PHTM    = vars[n_v+n_g+4] + vars_SS[n_v+n_g+4] # consumption of poor hand-to-mouth

        aggZ        = vars[n_v+n_g+n_p+1] + vars_SS[n_v+n_g+n_p+1] # aggregate Z
        V_Dot       = vars[nVars+1:nVars+n_v]
        g_Dot       = vars[nVars+n_v+1:nVars+n_v+n_g]
        aggZ_Dot    = vars[nVars+n_v+n_g+n_p+1]
        VEErrors    = vars[2*nVars+1:2*nVars+n_v]
        aggZ_Shock  = vars[2*nVars+nEErrors+1]

        # Prices
        w   = (1 - aalpha) * (K ^ aalpha) * n_SS ^ (-aalpha) *
            ((permanent == 0) ? exp(aggZ) ^ (1-aalpha) : 1.)
        r_a = aalpha * (K ^ (aalpha - 1)) * (((permanent == 0) ? exp(aggZ) : 1.) *
                                             n_SS) ^ (1 - aalpha) - ddelta
        # Auxiliary variables
        r_b_borr = r_b .+ borrwedge_SS

        # Other necessary objects
        r_a_grid     = repeat([r_a],I,J,N)

        l_grid       = permutedims(repeat(ones(N,1),1,I,J), [2 3 1])
        y_shock      = y .* exp.(kappa * aggZ * (y .- y_mean) ./ std(y))
        y_shock_mean = dot(y_shock, y_dist) #y_shock * y_dist
        y_shock      = y_shock ./ y_shock_mean .* y_mean
        y_grid       = permutedims(repeat(y_shock',1,I,J), [2 3 1])

        # ripped out
        @time c, s, u, d, d_g, s_g, c_g = eqcond_helper(V, I, J, I_g, J_g, N, chi0, chi1, chi2,
                                                        a_lb, ggamma, permanent, interp_decision,
                                                        ddeath, pam, aggZ, xxi, tau_I, w, l_grid,
                                                        y_grid, b_grid, r_b_grid, trans_grid,
                                                        a_grid,
                                                        daf_grid, dab_grid, dbf_grid, dbb_grid,
                                                        IcF, IcB, Ic0, IcFB, IcBF, IcBB, Ic00)

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
        perm_mult = (permanent == 0) ? rrho + ddeath : rrho + ddeath - (1 - ggamma) * aggZ
        hjbResidual = reshape(u, I*J*N, 1) + A * reshape(V, I*J*N, 1) + V_Dot + VEErrors -
            perm_mult * reshape(V,I*J*N,1)

        # KFE
        gResidual = g_Dot - gIntermediate[1:n_g, 1]

        K_out = 0.0
        if K_liquid == 1
            K_out = sum((vec(a_g_grid) .+ vec(b_g_grid)) .* gg .* vec(dab_g))
        else
            K_out = sum(vec(a_g_grid) .* gg .* vec(dab_g))
        end

        K_Residual = K_out - K

        r_b_out = 0.0
        if r_b_fix      == 1
            r_b_out          = r_b_SS
            r_b_Residual = r_b_out - r_b
        elseif r_b_phi  == 1
            r_b_out          = sum(vec(b_g_grid) .* gg .* vec(dab_g))
            r_b_Residual = r_b_out - B_SS * exp(1/pphi * (r_b - r_b_SS))
        elseif B_fix    == 1
            r_b_out          = 0.0
            r_b_Residual = r_b_out - b_save
        elseif K_liquid == 1
            r_b_out          = r_a_out - illiquid_wedge
            r_b_Residual = r_b_out - r_b
        end

        WHTM_indicator      = zeros(I_g,J_g,N)
        WHTM_indicator[b_g_0pos:b_g_0pos+1,a_g_0pos+2:end,:] .= 1.
        WHTM_out            = sum(vec(WHTM_indicator) .* gg .* vec(dab_g))
        C_WHTM_out          = sum(vec(WHTM_indicator) .* vec(c_g) .* gg .* vec(dab_g))

        PHTM_indicator      = zeros(I_g,J_g,N)
        PHTM_indicator[b_g_0pos:b_g_0pos+1,a_g_0pos:a_g_0pos+2:end,:] .= 1.
        PHTM_out            = sum(vec(PHTM_indicator) .* gg .* vec(dab_g))
        C_PHTM_out          = sum(vec(PHTM_indicator) .* vec(c_g) .* gg .* vec(dab_g))

        #NHTM_indicator      = zeros(I_g,J_g,N)
        #if permanent == 0
        #    NHTM_indicator[b_g_0pos+3:end,2:end,:] .= 1.
        #elseif permanent == 1
        #    NHTM_indicator[b_g_0pos+2:end,:,:] .= 1.
        #end
        #NHTM_out          = sum(vec(NHTM_indicator) .* gg .* vec(dab_g))
        #C_NHTM_out        = sum(vec(NHTM_indicator) .* vec(c_g) .* gg .* vec(dab_g)) / NHTM_out

        C_WHTM_Residual     = C_WHTM_out - C_WHTM
        C_PHTM_Residual     = C_PHTM_out - C_PHTM

        # Law of motion for aggregate tfp shock
        aggZ_Residual = aggZ_Dot - (-nnu_aggZ * aggZ + ssigma_aggZ * aggZ_Shock)

        # Return equilibrium conditions
        return [hjbResidual; gResidual; K_Residual; r_b_Residual; C_WHTM_Residual;
                C_PHTM_Residual; aggZ_Residual]
    end
=#
    l_grid       = permutedims(repeat(ones(N,1),1,I,J), [2 3 1])
    @inline function get_residuals(vars::Vector{T}) where {T<:Real}
        # Unpack variables
        V   = reshape(vars[1:n_v] .+ vars_SS[1:n_v], I, J, N)  # value function
        g   = vars[n_v + 1 : n_v + n_g] .+ vars_SS[n_v + 1 : n_v + n_g]    # distribution
        K   = vars[n_v + n_g + 1] + vars_SS[n_v + n_g + 1]    # aggregate capital
        r_b = vars[n_v + n_g + 2] + vars_SS[n_v + n_g + 2]

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

        VEErrors   = vars[2*nVars + 1 : 2 * nVars + n_v]
        aggZ_Shock = vars[2*nVars + nEErrors + 1]

        g_end     = (1 - sum(g .* dab_g_aux[1:end-1])) / dab_g[I_g,J_g,N]
        gg        = vcat(g, g_end)

        # Prices
        w   = (1 - aalpha) * (K ^ aalpha) * n_SS ^ (-aalpha) *
            ((permanent == 0) ? exp(aggZ) ^ (1-aalpha) : 1.)
        r_a = aalpha * (K ^ (aalpha - 1)) * (((permanent == 0) ? exp(aggZ) : 1.) *
                                             n_SS) ^ (1 - aalpha) - ddelta
        # Auxiliary variables
        r_b_borr = r_b .+ borrwedge_SS

        # Other necessary objects
        r_a_grid     = repeat([r_a],I,J,N)

        y_shock      = y .* exp.(kappa * aggZ * (y .- y_mean) ./ std(y))
        y_shock_mean = dot(y_shock, y_dist) #y_shock * y_dist
        y_shock      = y_shock ./ y_shock_mean .* y_mean
        y_grid       = permutedims(repeat(y_shock',1,I,J), [2 3 1])

        # ripped out
        @show "eqcond_helper"
        @time c, s, u, d, d_g, s_g, c_g = eqcond_helper(V, I, J, I_g, J_g, N, chi0, chi1, chi2,
                                                        a_lb, ggamma, permanent, interp_decision,
                                                        ddeath, pam, aggZ, xxi, tau_I, w, l_grid,
                                                        y_grid, b_grid, r_b_grid, trans_grid,
                                                        a_grid,
                                                        daf_grid, dab_grid, dbf_grid, dbb_grid,
                                                        IcF, IcB, Ic0, IcFB, IcBF, IcBB, Ic00)

        # Derive transition matrices
        @show "transition_deriva"
        @time aa, bb, aau, bbu = transition_deriva(I_g, J_g, N, I, J, permanent, ddeath, pam,
                                                   xxi, w, chi0, chi1, chi2, a_lb, l_grid,
                                                   l_g_grid, y_grid, y_g_grid, d,
                                             dab_grid, daf_grid, dab_g_grid, daf_g_grid, dbb_grid,
                                             dbf_grid, dbb_g_grid, dbf_g_grid, d_g, a_grid,
                                             a_g_grid, s, s_g, r_a_grid, r_a_g_grid, aggZ)

        # full transition matrix
        A = aa + bb + cc

        #----------------------------------------------------------------
        # KFE
        #----------------------------------------------------------------
        AT = (aau + bbu + ccu)'
        #AT = AT'
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
        a_g_grid_aux = reshape(a_g_grid, I_g*J_g*N, 1)
        a_save       = sum(gIntermediate .* dab_g_aux .* a_g_grid_aux)

        b_g_grid_aux = reshape(b_g_grid, I_g*J_g*N, 1)
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
        perm_mult = (permanent == 0) ? rrho + ddeath : rrho + ddeath - (1 - ggamma) * aggZ
        hjbResidual = reshape(u, I*J*N, 1) + A * reshape(V, I*J*N, 1) + V_Dot + VEErrors -
            perm_mult * reshape(V,I*J*N,1)

        # KFE
        gResidual = g_Dot - gIntermediate[1:n_g, 1]

        K_out = 0.0
        if K_liquid == 1
            K_out = sum((vec(a_g_grid) .+ vec(b_g_grid)) .* gg .* vec(dab_g))
        else
            K_out = sum(vec(a_g_grid) .* gg .* vec(dab_g))
        end

        K_Residual = K_out - K
        r_b_out = 0.0
        r_b_Residual = 0.0
        if r_b_fix      == 1
            r_b_out      = r_b_SS
            r_b_Residual = r_b_out - r_b
        elseif r_b_phi  == 1
            r_b_out      = sum(vec(b_g_grid) .* gg .* vec(dab_g))
            r_b_Residual = r_b_out - B_SS * exp(1/pphi * (r_b - r_b_SS))
        elseif B_fix    == 1
            r_b_out      = 0.0
            r_b_Residual = r_b_out - b_save
        elseif K_liquid == 1
            r_b_out      = r_a_out - illiquid_wedge
            r_b_Residual = r_b_out - r_b
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

            PHTM_indicator      = zeros(I_g,J_g,N)
            PHTM_indicator[b_g_0pos:b_g_0pos+1,a_g_0pos:a_g_0pos+2:end,:] .= 1.
            PHTM_out            = sum(vec(PHTM_indicator) .* gg .* vec(dab_g))
            C_PHTM_out          = sum(vec(PHTM_indicator) .* vec(c_g) .* gg .* vec(dab_g))

            if permanent == 0
                NHTM_indicator      = zeros(I_g,J_g,N)
                NHTM_indicator[b_g_0pos+2:end,a_g_0pos+2:end,:] .= 1.

                NHTM_indicator      = zeros(I_g,J_g,N)
                NHTM_indicator[b_g_0pos+3:end,2:end,:] .= 1.
                NHTM_out            = sum(vec(NHTM_indicator) .* gg .* vec(dab_g))
                C_NHTM_out          = sum(vec(NHTM_indicator) .* vec(c_g) .* gg .*
                                          vec(dab_g))/NHTM_out

            elseif permanent == 1
                NHTM_indicator      = zeros(I_g,J_g,N)
                NHTM_indicator[b_g_0pos+2:end,:,:] .= 1.
                NHTM_out            = sum(vec(NHTM_indicator) .* gg .* vec(dab_g))
                C_NHTM_out          = sum(vec(NHTM_indicator) .* vec(c_g) .* gg .*
                                          vec(dab_g))/NHTM_out
            end
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

    # equilibrium conditions
    #=
    if aggregate_variables == 1
        println("Agg")
        out = get_residuals_agg(zeros(Float64, 2 * nVars + nEErrors + 1))
        @time get_residuals_agg(zeros(Float64, 2 * nVars + nEErrors + 1))
    elseif distributional_variables == 1
        println("Dist")
        out = get_residuals_dist(zeros(Float64, 2 * nVars + nEErrors + 1))
        @time get_residuals_dist(zeros(Float64, 2 * nVars + nEErrors + 1))
    elseif distributional_variables_1 == 1
        println("Dist_1")
        out = get_residuals_dist1(zeros(Float64, 2 * nVars + nEErrors + 1))
        @time get_residuals_dist1(zeros(Float64, 2 * nVars + nEErrors + 1))
    else
    #   vResidual = [hjbResidual; gResidual; K_Residual; r_b_Residual; aggZ_Residual]
    end
    =#

    out = get_residuals(zeros(Float64, 2 * nVars + nEErrors + 1))
    #JLD2.jldopen("eqcond_before_factoring.jld2", true, true, true, IOStream) do file
    #    file["residuals"] = out
    #end
    test_out = load("eqcond_before_factoring.jld2", "residuals")
    @assert test_out == out

    @time get_residuals(zeros(Float64, 2 * nVars + nEErrors + 1))
    error()
    derivs = ForwardDiff.jacobian(get_residuals, zeros(Float64, 2 * nVars + nEErrors + 1))

    #####
    #nstates = n_states(m)
    #n_s_exp = n_shocks_expectational(m)
    #n_s_exo = n_shocks_exogenous(m)
    #vars = zeros(Float64, 2 * nstates + n_s_exp + n_s_exo)
    ####

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
