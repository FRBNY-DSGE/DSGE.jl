using SparseArrays
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
    pam    = m[:pam].value
    xxi    = m[:xxi].value
    ggamma = m[:ggamma].value
    tau_I  = m[:tau_I].value
    trans  = m[:trans].value

    # Set liquid rates
    r_b      = m[:r_b_SS].value
    r_b_borr = m[:r_b_borr_SS].value
    lambda   = get_setting(m, :lambda)
    K_liquid = get_setting(m, :K_liquid)

    # Output flags
    aggregate_variables        = get_setting(m, :aggregate_variables)
    distributional_variables   = get_setting(m, :distributional_variables)
    distributional_variables_1 = get_setting(m, :distributional_variables_1)

    # Approximation parameters
    maxit_HJB = get_setting(m, :maxit_HJB)
    crit_HJB  = get_setting(m, :crit_HJB)
    Delta     = get_setting(m, :Delta)

    maxit_HIS = get_setting(m, :maxit_HIS)
    crit_HIS  = get_setting(m, :crit_HIS)
    start_HIS = get_setting(m, :start_HIS)

    maxit_KFE = get_setting(m, :maxit_KFE)
    crit_KFE  = get_setting(m, :crit_KFE)
    Δ_KFE     = get_setting(m, :Delta_KFE)

    maxit_KL  = get_setting(m, :maxit_KL)
    crit_KL   = get_setting(m, :crit_KL)
    relax_KL  = get_setting(m, :relax_KL)

    I      = get_setting(m, :I)
    J      = get_setting(m, :J)
    I_g    = get_setting(m, :I_g)
    J_g    = get_setting(m, :J_g)
    N      = get_setting(m, :N)

    # Read in grids
    a      = get_setting(m, :a)
    b      = get_setting(m, :b)
    a_g    = get_setting(m, :a_g)
    b_g    = get_setting(m, :b_g)
    y      = real(vec(get_setting(m, :y)))
    y_dist = get_setting(m, :y_dist)
    y_mean = get_setting(m, :y_mean)
    KL     = get_setting(m, :KL_0)
    interp_decision = get_setting(m, :interp_decision)

    # Compute prices associated with initial guess of KL
    w	= (1 - aalpha) * (KL ^ aalpha)
    r_a	= aalpha * (KL ^ (aalpha - 1)) - ddelta

    # Construct problem functions
    util, deposit, cost = construct_problem_functions(ggamma, chi0, chi1, chi2, a_lb)

    # Initialize difference grids
    a_g_grid, b_g_grid, dab_g_tilde_grid, dab_g_tilde = set_grids(a, b, a_g, b_g, N)

    dab_g_tilde_vec = repeat(dab_g_tilde, N)
    r_b_vec         = r_b .* (b .>= 0) + r_b_borr .* (b .< 0)

    b_vec = repeat(repeat(b, outer=J), outer=N)
    a_vec = repeat(repeat(a, inner=I), outer=N)
    y_vec = repeat(y, inner=I*J)

    y_grid   = reshape(repeat(y, inner=I*J), I, J, N)
    y_g_grid = reshape(repeat(y, inner=I_g*J_g), I_g, J_g, N)

    # Initialization for consumption and value function
    c_0 = (1-xxi) * w * y_vec .+ (r_a + ddeath*pam) .* a_vec +
        (r_b_borr + ddeath * pam) .* b_vec .+ trans - tau_I * w * y_vec

    V_0	= reshape(1/(rrho + ddeath) * util.(c_0), I, J, N)

    # Initial distribution
    gg0 = zeros(I_g, J_g, N)
    gg0[b .== 0, 1, :] = vec(y_dist)
    gg0 = vec(gg0)
    gg0 = (gg0 ./ sum(gg0)) ./ dab_g_tilde_vec # ensure integration to 1
    gg  = gg0

    AT  = Array{Float64}(undef, 0, 0)
    Vn1 = Array{Float64}(undef, I, J, N)
    Vn  = Array{Float64}(undef, I, J, N)

    gg_tilde = Array{Float64}(undef, 0, 0)

    g   = Array{Float64}(undef, 0, 0)
    c   = Array{Float64,3}(undef, I, J, N)
    s   = Array{Float64,3}(undef, I, J, N)
    u   = Array{Float64,3}(undef, I, J, N)
    d   = Array{Float64,3}(undef, I, J, N)
    c_g = Array{Float64,3}(undef, I_g, J_g, N)

    K_supply, L_supply, KL_supply = 0.0, 0.0, 0.0

    #----------------------------------------------------------------
    # Iterate on KL to find steady state
    #----------------------------------------------------------------
    for ii = 1 : maxit_KL

        # Derive aggregates, given KL
	    w   = (1 - aalpha) * (KL ^ aalpha)
 	    r_a = aalpha * (KL ^ (aalpha - 1)) - ddelta

	    # Store current value function
	    Vn = V_0

	    #----------------------------------------------------------------
	    # Solve HJB
	    #----------------------------------------------------------------
        @time for nn = 1 : maxit_HJB
            c, s, d = solve_hjb(Vn, a_lb, ggamma, ddeath, pam, trans, xxi,
                                tau_I, 0.0, w, r_b_vec, y, a, b, cost, util, deposit)
            u = util.(c)

            # Interpolate
            d_g = reshape(interp_decision * vec(d), I_g, J_g, N)
            s_g = reshape(interp_decision * vec(s), I_g, J_g, N)
            c_g = reshape(interp_decision * vec(c), I_g, J_g, N)

            A, AT = transition(ddeath, pam, xxi, w, a_lb, 0.0, d, d_g, s, s_g,
                               r_a, a, a_g, b, b_g, y, cost)

            #------------------------------------------------------------
            # Update value function, Howard improvement steps
            #------------------------------------------------------------
            Vn1 = update_value_fn(A, Vn, lambda, u, Delta, rrho, ddeath)

            if nn >= start_HIS
                for jj = 1:maxit_HIS
                    Vn2    = update_value_fn(A, Vn1, lambda, u, Delta, rrho, ddeath)
                    VHIS_Δ = Vn2 - Vn1
                    Vn1    = Vn2
                    dist   = maximum(abs.(VHIS_Δ))
                    if dist < crit_HIS
                        break
                    end
                end
            end

            # Check for convergence
            V_Δ  = Vn1 - Vn
            Vn   = Vn1
            dist = maximum(abs.(V_Δ))
            @show dist

            if dist < crit_HJB
                 println("Value Function Converged, Iteration = ", nn)
                break
            end
        end

        # Store final value function
        V_0 = Vn

        #----------------------------------------------------------------
        # Compute stationary distribution and update aggregate KL
        #----------------------------------------------------------------

        # Find new stationary distribution associated with decision rules
        A   = AT
        λ0  = lambda - diagm(0 => diag(lambda)) # transition matrix with diagonal killed
        λ0p = λ0'

        gg_tilde = dab_g_tilde_vec .* gg
        gg1      = Array{Float64,2}(undef, I_g * J_g, N)
        g        = Array{Float64,3}(undef, I_g,  J_g, N)

        K_supply, L_supply, KL_supply = 0.0, 0.0, 0.0

        #----------------------------------------------------------------
        # Iteration to solve Kolmogorov Forward equation
        #----------------------------------------------------------------
        @time for nn = 1:maxit_KFE

            gg_tilde  = dab_g_tilde_vec .* gg
            gg_tilde2 = reshape(gg_tilde, I_g*J_g, N)
            gg1       = Array{Float64,2}(undef, I_g * J_g, N)

            for kk = 1:N
                Ak = A[1+(kk-1)*(I_g*J_g):kk*(I_g*J_g), 1+(kk-1)*(I_g*J_g):kk*(I_g*J_g)]

                death_inflow = zeros(I_g, J_g)
                death_inflow[b .== 0, 1] .= sum(gg_tilde2[:,kk])

                gk_sum = vec(sum(broadcast(*, λ0p[kk,:], gg_tilde2'), dims=1))

                gg1[:,kk] = add_diag(add_diag(-Δ_KFE*Ak, 1.0),
                                     -Δ_KFE * (lambda[kk,kk] - ddeath))' \
                                         (gg_tilde2[:,kk] + Δ_KFE * gk_sum +
                                          Δ_KFE * ddeath * vec(death_inflow))
            end

            gg11 = (vec(gg1) ./ sum(gg1)) ./ dab_g_tilde_vec
            dist = maximum(abs.(gg11-gg))

            if mod(nn,100) == 0
                println("Iteration ", nn, ", distance is ", dist)
            end
            if dist < crit_KFE
                gg = gg11
                g  = reshape(gg, I_g, J_g, N)

                println("Distribution Converged, Iteration = ", nn)
                break
            end
            gg = gg11
            g  = reshape(gg, I_g, J_g, N)
        end

        #------------------------------------------------------------
        # Update guess of Capital, Labor Supply
        #------------------------------------------------------------
        ggdab_g = g .* dab_g_tilde_grid

        K_supply  = sum((a_g_grid .+ (K_liquid ? b_g_grid : 0.0)) .* ggdab_g)
        L_supply  = sum(y_g_grid .* ggdab_g)
        KL_supply = K_supply / L_supply

        # Check for convergence and update
        KL_Δ = KL - KL_supply
        println("The current gap is ", KL_Δ)

        if abs(KL_Δ) > crit_KL
            KL = relax_KL * KL_supply + (1 - relax_KL) * KL
        else
            println("I have found the steady state, Iteration = ", ii)
            break
        end
    end

    ccu                  = kron(lambda, my_speye(I_g * J_g))
    A                    = AT + ccu
    dab_small            = reshape(dab_g_tilde_grid[:,:,1], I_g*J_g, 1)
    loc                  = findall(!iszero, b .== 0)
    dab_small            = dab_small ./ dab_small[loc] * ddeath
    dab_small[loc]      .= 0.0
    death_process        = -ddeath * my_speye(I_g*J_g)
    death_process[loc,:] = vec(dab_small)
    death_process        = kron(my_speye(N), death_process)

    g_new   = A' * gg_tilde
    g_new   = (1.0 ./ dab_g_tilde_vec) .* g_new + death_process * gg
    g_new   = reshape(g_new ,I_g, J_g, N)
    g_new_a = sum(sum(g_new, dims=1), dims=3)
    #compute_savings()

    # Rename variables in steady state
    m[:KL_SS]           = real(KL)
    m[:r_a_SS]          = aalpha * (m[:KL_SS].value ^ (aalpha - 1)) - ddelta
    m[:w_SS]            = (1 - aalpha) * (m[:KL_SS].value ^ aalpha)
    m[:K_SS]            = real(m[:KL_SS].value * L_supply)
    m[:u_SS]            = u
    m[:c_SS]            = c
    m[:d_SS]            = d
    m[:V_SS]            = Vn
    m[:g_SS]            = g
    m[:C_SS]            = sum(g .* c_g .* dab_g_tilde_grid)
    m[:C_Var_SS]        = sum(g .* log.(c_g).^2 .* dab_g_tilde_grid) -
                                 sum(g .* log.(c_g) .* dab_g_tilde_grid)^2
    m[:I_SS]            = dot(g_new_a, a_g)
    m[:B_SS]            = sum(g .* b_g_grid .* dab_g_tilde_grid)
    m[:Y_SS]            = real((K_supply ^ aalpha) * (L_supply ^ (1 - aalpha)))
    m[:n_SS]            = real(L_supply)

    r_b_g_grid          = r_b .* (b_g_grid .>= 0) + r_b_borr .* (b_g_grid .< 0)
    m[:earn_SS]         = log.((1-tau_I) * w * y_g_grid + b_g_grid .*
                               (r_b_g_grid .+ ddeath*pam) .+ trans .+ a_g_grid .*
                               (r_a .+ ddeath*pam))

    m[:earn_Var_SS]     = sum(g .* m[:earn_SS].value .^ 2 .* dab_g_tilde_grid) -
                                sum(g .* m[:earn_SS].value .* dab_g_tilde_grid)^2

    a_g_0pos = get_setting(m, :a_g_0pos)
    b_g_0pos = get_setting(m, :b_g_0pos)

    ###
    # Consumption by hand-to-mouth status
    ###
    WHTM_indicator      = zeros(I_g, J_g, N)
    WHTM_indicator[b_g_0pos:b_g_0pos+1, a_g_0pos+2:end,:] .= 1.0

    m[:WHTM_indicator] = WHTM_indicator
    m[:WHTM_SS]        = sum(g[:] .* WHTM_indicator[:] .* dab_g_tilde_grid[:])
    m[:C_WHTM_SS]      = sum(WHTM_indicator[:] .* c_g[:] .* g[:] .* dab_g_tilde_grid[:])

    PHTM_indicator     = zeros(I_g, J_g, N)
    PHTM_indicator[b_g_0pos:b_g_0pos+1,a_g_0pos:a_g_0pos+2:end,:] .= 1.0
    m[:PHTM_indicator] = PHTM_indicator
    m[:PHTM_SS]        = sum(g[:] .* PHTM_indicator[:] .* dab_g_tilde_grid[:])
    m[:C_PHTM_SS]      = sum(c_g[:] .* g[:] .* PHTM_indicator[:] .* dab_g_tilde_grid[:])

    NHTM_indicator      = zeros(I_g,J_g,N)
    NHTM_indicator[(b_g_0pos+3):end,2:end,:] .= 1.0

    m[:NHTM_indicator] = NHTM_indicator
    m[:NHTM_SS]        = sum(g[:] .* NHTM_indicator[:] .* dab_g_tilde_grid[:])
    m[:C_NHTM_SS]      = sum(c_g[:] .* g[:] .* NHTM_indicator[:] .* dab_g_tilde_grid[:]) /
        m[:NHTM_SS].value

    n_v = get_setting(m, :n_v)
    n_g = get_setting(m, :n_g)
    n_p = get_setting(m, :n_p)
    n_Z = get_setting(m, :n_Z)

    # Collect variables into vector
    vars_SS                  = zeros(n_v + n_g + n_p + n_Z)
    vars_SS[1:n_v]           = vec(m[:V_SS].value)
    gg_SS                    = vec(m[:g_SS].value)
    m[:gg_SS]                = gg_SS

    vars_SS[n_v+1:n_v+n_g] = gg_SS[1:I_g*J_g*N-1]
    vars_SS[n_v+n_g+1]     = m[:K_SS].value
    vars_SS[n_v+n_g+2]     = m[:r_b_SS].value

    if aggregate_variables == 1
        vars_SS[n_v+n_g+3] = m[:Y_SS].value
        vars_SS[n_v+n_g+4] = m[:C_SS].value
    elseif distributional_variables == 1
        vars_SS[n_v+n_g+3] = m[:C_Var_SS].value
        vars_SS[n_v+n_g+4] = m[:earn_Var_SS].value
    elseif distributional_variables_1 == 1
        vars_SS[n_v+n_g+3] = m[:C_WHTM_SS].value
        vars_SS[n_v+n_g+4] = m[:C_PHTM_SS].value
    end
    vars_SS[n_v+n_g+n_p+1,1] = 0.0

    # Compute illiquid wedge
    m[:illiquid_wedge] = m[:r_a_SS].value - m[:r_b_SS].value

    # Save SS variables
    m[:vars_SS] = vars_SS

    return m
end
