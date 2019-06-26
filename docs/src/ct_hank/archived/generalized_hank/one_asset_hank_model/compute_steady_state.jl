## Solves the (non-stochastic) steady state of one-asset HANK model

function compute_steady_state(grids::Dict{Symbol, Any},
                              params::Dict{Symbol, Float64},
                              init_params::Dict{Symbol, Float64},
                              approx_params::Dict{Symbol, Any})

     # Read in parameters
    coefrra = params[:coefrra]
    frisch  = params[:frisch]
    meanlabeff = params[:meanlabeff]
    maxhours = params[:maxhours]

    ceselast = params[:ceselast]
    labtax = params[:labtax]                        # marginal tax rate
    govbondtarget = params[:govbondtarget]          # multiple of quarterly GDP
    labdisutil = params[:labdisutil]
    lumptransferpc = params[:lumptransferpc]        # 6% of quarterly GDP in steady state

    IterateR = approx_params[:IterateR]
    IterateRho = approx_params[:IterateRho]


    # Read in grids
    a   = grids[:a]
    g_z = grids[:g_z]
    zz  = grids[:zz]
    ymarkov_combined = grids[:ymarkov_combined]

    I   = length(a)
    J   = length(g_z)
    amax = maximum(a)
    amin = minimum(a)

    # Read in initial rates
    r = init_params[:r0]
    rmin = init_params[:rmin]
    rmax = init_params[:rmax]

    rrho = init_params[:rrho0]
    rhomin = init_params[:rhomin]
    rhomax = init_params[:rhomax]

    # Read in approximation parameters
    Ir = approx_params[:Ir]
    maxit_HJB = approx_params[:maxit_HJB]
    tol_HJB = approx_params[:tol_HJB]
    Delta_HJB = approx_params[:Delta_HJB]
    maxit_KFE = approx_params[:maxit_KFE]
    tol_KFE = approx_params[:tol_KFE]
    Delta_KFE = approx_params[:Delta_KFE]
    niter_hours = approx_params[:niter_hours]
    crit_S = approx_params[:crit_S]

    # Initializing equilibrium objects
    m_SS = (ceselast - 1)/ceselast
    w = w_SS = m_SS
    N_SS, Y_SS, B_SS, profit_SS, profshare, lumptransfer =
    calculate_SS_equil_vars(zz, m_SS, meanlabeff, lumptransferpc, govbondtarget)

    # Initialize matrices for finite differences
    Vaf = zeros(Complex128, I, J)
    Vab = zeros(Complex128, I, J)

    cf  = similar(Vaf)
    hf  = similar(Vaf)
    sf  = similar(Vaf)

    cb  = similar(Vab)
    hb  = similar(Vab)
    sb  = similar(Vab)

    c0 = similar(cb)

    A = zeros(Complex128, 0)
    Aswitch = kron(ymarkov_combined, speye(Complex128, I))

    # Initialize steady state variables
    vars_SS = OrderedDict{Symbol, Any}()
    V   = similar(Vaf)
    u   = similar(Vaf)
    s   = similar(Vaf)

    util, income, labor = construct_household_problem_functions(V, w, params)

    # Setting up forward/backward difference grids for a
    dazf, dazb, azdelta, aa, adelta, azdelta, azdelta_mat = initialize_diff_grids(zz, a)

    for ir = 1:Ir
        c  = zeros(Complex128, I, J)
        h  = fill(complex(1/3), I, J)
        h0 = ones(Complex128, I, J)

        # initial guess
        v = similar(h)
        for i in eachindex(h)
            inc  = income(h[i], zz[i], profshare[i], lumptransfer, r, aa[i])
            v[i] = util(inc, h[i])/rrho
        end

        # iterate HJB
        for ihjb = 1:maxit_HJB
            V = v

            Vaf, Vab, cf, hf, cb, hb = construct_initial_diff_matrices(V, Vaf, Vab, income,
                                                                       labor, h, h0, zz,
                                                                       profshare,
                                                                       lumptransfer, amax, amin,
                                                                       coefrra, r, dazf,
                                                                       dazb, maxhours)

            cf, hf, cb, hb, c0, h0 = hours_iteration(income, labor, zz, profshare,
                                                     lumptransfer, aa,
                                                     coefrra, r, cf, hf, cb, hb, c0, h0,
                                                     maxhours, niter_hours)

            for i in eachindex(h0)
                c0[i] = income(h0[i], zz[i], profshare[i], lumptransfer, r, aa[i])

                sf[i] = income(hf[i], zz[i], profshare[i], lumptransfer, r, aa[i]) - cf[i]
                sb[i] = income(hb[i], zz[i], profshare[i], lumptransfer, r, aa[i]) - cb[i]
            end

            Vaf[I, :] = cf[I, :] .^ (-coefrra)
            Vab[1, :] = cb[1, :] .^ (-coefrra)
            Va0 = (c0).^(-coefrra)

            Vf = similar(cf); Vb = similar(cb); V0 = similar(c0)

            for i in eachindex(cf)
                Vf[i] = util(cf[i], hf[i])
                Vb[i] = util(cb[i], hb[i])
                V0[i] = util(c0[i], h0[i])
            end

            # Compute upwinded value function difference, and indicator matrices
            dV_upwind, If, Ib, I0 = upwind_value_function(Vaf, Vab, Va0, Vf, Vb, V0, cf, cb, c0, sf, sb)

            # Compute labor, consumption, savings, and flow utility based on upwind scheme
            h = hf .* If + hb .* Ib + h0 .* I0
            c = cf .* If + cb .* Ib + c0 .* I0
            s = sf .* If + sb .* Ib
            u = similar(c)
            for i in eachindex(c)
                u[i] = util(c[i], h[i])
            end

            A = upwind_matrix(Aswitch, sf, sb, dazf, dazb, I, 2; exact = false, If = If, Ib = Ib, I0 = I0)
            V_stacked = solve_hjb(A, rrho, Delta_HJB, vec(u), vec(V))
            V = reshape(V_stacked, I, 2)

            Vchange = V - v
            v = V

            err_HJB = maximum(abs.(Vchange))
            if err_HJB < tol_HJB
                break
            end
        end

        g = solve_KFE(A, a, g_z, azdelta, azdelta_mat, maxit_KFE = maxit_KFE, tol_KFE = tol_KFE,
                      Delta_KFE = Delta_KFE)

        @assert false

        N_SS, Y_SS, B_SS, profit_SS, profshare, lumptransfer, bond_err =
        calculate_SS_equil_vars(zz, h, g, azdelta, aa, m_SS, meanlabeff,
                                lumptransferpc, govbondtarget)


        r, rmin, rmax, rrho, rhomin, rhomax, clear_cond =
        check_bond_market_clearing(bond_err, crit_S, r, rmin, rmax, rrho, rhomin, rhomax,
                                   IterateR, IterateRho)

        if clear_cond
            vars_SS[:V_SS] = real(V)
            vars_SS[:inflation_SS] = 0.
            vars_SS[:g_SS] = real(g)
            vars_SS[:r_SS] = r
            vars_SS[:u_SS] = real(u)
            vars_SS[:c_SS] = real(c)
            vars_SS[:h_SS] = real(h)
            vars_SS[:s_SS] = real(s)
            vars_SS[:rnom_SS] = vars_SS[:r_SS] + vars_SS[:inflation_SS]
            vars_SS[:B_SS] = sum(vec(vars_SS[:g_SS]) .* vec(aa) .* vec(azdelta))
            vars_SS[:N_SS] = real(sum(vec(zz) .* vec(vars_SS[:h_SS]) .*
                                      vec(vars_SS[:g_SS]) .* vec(azdelta)))
            vars_SS[:Y_SS] = vars_SS[:N_SS]
            vars_SS[:m_SS] = (ceselast - 1)/ceselast
            vars_SS[:w_SS] = vars_SS[:m_SS]
            vars_SS[:profit_SS] = (1 - vars_SS[:m_SS]) .* vars_SS[:Y_SS]
            vars_SS[:C_SS] = sum(vec(vars_SS[:c_SS]) .* vec(vars_SS[:g_SS]) .*
                                 vec(azdelta))
            vars_SS[:T_SS] = real(lumptransfer)
            vars_SS[:G_SS] = labtax * vars_SS[:w_SS] * vars_SS[:N_SS] - vars_SS[:T_SS] -
                             vars_SS[:r_SS] * vars_SS[:B_SS]
            vars_SS[:rrho] = rrho
            break
        end

    end

    return vars_SS
end
