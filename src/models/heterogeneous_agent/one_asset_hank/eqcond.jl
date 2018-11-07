"""
``
eqcond(m::OneAssetHANK)
```

Expresses the equilibrium conditions in canonical form using Γ0, Γ1, C, Ψ, and Π matrices.

### Outputs

* `Γ0` (`n_states` x `n_states`) holds coefficients of current time states.
* `Γ1` (`n_states` x `n_states`) holds coefficients of lagged states.
* `C`  (`n_states` x `1`) is a vector of constants
* `Ψ`  (`n_states` x `n_shocks_exogenous`) holds coefficients of iid shocks.
* `Π`  (`n_states` x `n_states_expectational`) holds coefficients of expectational states.
"""
function eqcond(m::OneAssetHANK)
    nstates = n_states(m)
    n_s_exp = n_shocks_expectational(m)
    n_s_exo = n_shocks_exogenous(m)

    x       = zeros(Float64, 2 * nstates + n_s_exp + n_s_exo)

    # Read in parameters, settings
    niter_hours       = get_setting(m, :niter_hours)
    n_v               = get_setting(m, :n_jump_vars)
    n_g               = get_setting(m, :n_state_vars)
    ymarkov_combined  = get_setting(m, :ymarkov_combined)
    zz                = get_setting(m, :zz)
    a                 = get_setting(m, :a)
    I, J              = size(zz)
    amax              = get_setting(m, :amax)
    amin              = get_setting(m, :amin)
    aa                = repeat(a, 1, J)
    A_switch          = kron(ymarkov_combined, speye(Float64, I))
    daf, dab, azdelta = initialize_diff_grids(a, I, J)

    # Set up steady state parameters
    V_ss = m[:V_ss].value
    g_ss = m[:g_ss].value
    G_ss = m[:G_ss].value
    w_ss = m[:w_ss].value
    N_ss = m[:N_ss].value
    C_ss = m[:C_ss].value
    Y_ss = m[:Y_ss].value
    B_ss = m[:B_ss].value
    r_ss = m[:r_ss].value
    ρ_ss = m[:ρ_ss].value
    θ_MP = m[:θ_MP].value
    σ_MP = m[:σ_MP].value

    h_ss              = reshape(m[:h_ss].value, I, J)
    ceselast          = m[:ceselast].value
    inflation_ss      = m[:inflation_ss].value
    maxhours          = m[:maxhours].value
    govbcrule_fixnomB = m[:govbcrule_fixnomB].value
    priceadjust       = m[:priceadjust].value
    taylor_inflation  = m[:taylor_inflation].value
    taylor_outputgap  = m[:taylor_outputgap].value
    meanlabeff        = m[:meanlabeff].value

    # Necessary for construction of household problem functions
    labtax            = m[:labtax].value
    coefrra           = m[:coefrra].value
    frisch            = m[:frisch].value
    labdisutil        = m[:labdisutil].value

    # Note: make TFP a setting?
    TFP               = 1.0

    @inline function get_residuals(x::Vector{T}) where {T<:Real}
        # Prepare steady state deviations
        V           = reshape(x[1:n_v - 1]       + V_ss, I, J)
        inflation   = x[n_v]                     + inflation_ss
        gg          = x[n_v + 1 : n_v + n_g - 1] + g_ss[1:end-1]
        MP          = x[n_v + n_g]
        w           = x[n_v + n_g + 1]           + w_ss
        hours       = x[n_v + n_g + 2]           + N_ss
        consumption = x[n_v + n_g + 3]           + C_ss
        output      = x[n_v + n_g + 4]           + Y_ss
        assets      = x[n_v + n_g + 5]           + B_ss

        # Set up system of time derivatives and expectational errors, etc.
        V_dot           = x[nstates + 1 : nstates + n_v - 1]
        inflation_dot   = x[nstates + n_v]
        g_dot           = x[nstates + n_v + 1 : nstates + n_v + n_g - 1]
        mp_dot          = x[nstates + n_g + n_v]
        VEErrors        = x[2*nstates + 1 : 2*nstates + n_v - 1]
        inflation_error = x[2*nstates + n_v]
        mp_shock        = x[2*nstates + n_v + 1]

        g_end = (1 - dot(gg, azdelta[1:end-1])) ./ azdelta[end]
        g     = vcat(gg, g_end)

        #-----------------------------------------------------------------
        # Get equilibrium values, given steady state values
        #-----------------------------------------------------------------
        normalized_wage = w / TFP
        profshare       = (zz ./ meanlabeff) .* ((1. - normalized_wage) .* output)
        r_nominal       = r_ss + taylor_inflation * inflation + taylor_outputgap * (log(output)-log(Y_ss)) + MP
        r               = r_nominal - inflation
        lumptransfer    = labtax * w * hours - G_ss -
            (r_nominal - (1-govbcrule_fixnomB) * inflation) * assets

        #-----------------------------------------------------------------
        # Compute one iteration of the HJB
        #-----------------------------------------------------------------

        # Get flow utility, income, and labor hour functions
        util, income, labor = construct_household_problem_functions(V, w, coefrra, frisch, labtax, labdisutil)

        # Initialize other variables, using V to ensure everything is a dual number
        Vaf = copy(V)
        Vab = copy(V)
        cf  = similar(V)
        hf  = similar(V)
        sf  = similar(V)
        cb  = similar(V)
        hb  = similar(V)
        sb  = similar(V)
        c0  = similar(V)
        h0  = similar(V)
        for j=1:J, i=1:I
            h0[i,j] = h_ss[i,j]
        end

        #-----------------------------------------------------------------
        # Construct Initial Difference Matrices
        #-----------------------------------------------------------------
        for j=1:J, i=1:I
            if i==I
                Vaf[end,j] = income(h_ss[end,j], zz[end,j], profshare[end,j], lumptransfer, r, amax)^(-coefrra)
                Vab[i, j]  = (V[i,j] - V[i-1,j]) / dab[i]
            elseif i==1
                Vaf[i, j] = (V[i+1,j] - V[i,j]) / daf[i]
                Vab[1, j] = income(h_ss[1, j], zz[1, j], profshare[1, j], lumptransfer, r, amin)^(-coefrra)
            else
                Vaf[i, j] = (V[i+1,j] - V[i,j]) / daf[i]
                Vab[i, j] = (V[i,j] - V[i-1,j]) / dab[i]
            end

            cf[i,j] = Vaf[i,j] ^ (-1 / coefrra)
            cb[i,j] = Vab[i,j] ^ (-1 / coefrra)

            hf[i,j] = min(labor(zz[i,j], Vaf[i,j]), maxhours)
            hb[i,j] = min(labor(zz[i,j], Vab[i,j]), maxhours)
        end

        #-----------------------------------------------------------------
        # Hours Iteration
        #-----------------------------------------------------------------
        for ih=1:niter_hours
            for j=1:J, i=1:I
                if i == I
                    cf[end, j]  = income(hf[end, j], zz[end, j], profshare[end, j], lumptransfer, r, aa[end, j])
                    hf[end, j]  = min(labor(zz[end, j], cf[end, j]^(-coefrra)), maxhours)

                    if (ih == niter_hours) Vaf[end, j] = cf[I, j] ^ (-coefrra) end

                elseif i == 1
                    cb[1, j]  = income(hb[1, j], zz[1, j], profshare[1, j], lumptransfer, r, aa[1, j])
                    hb[1, j]  = min(labor(zz[1, j], cb[1, j]^(-coefrra)), maxhours)

                    if (ih == niter_hours) Vab[1, j] = cb[1, j] ^ (-coefrra) end
                end

                c0[i, j] = income(h0[i, j], zz[i, j], profshare[i, j], lumptransfer, r, aa[i, j])
                h0[i, j] = min(labor(zz[i, j], c0[i, j]^(-coefrra)), maxhours)
            end
        end

        for i in eachindex(h0)
            c0[i] = income(h0[i], zz[i], profshare[i], lumptransfer, r, aa[i])
            sf[i] = income(hf[i], zz[i], profshare[i], lumptransfer, r, aa[i]) - cf[i]
            sb[i] = income(hb[i], zz[i], profshare[i], lumptransfer, r, aa[i]) - cb[i]
        end

        A, u, h, c, s  = upwind(util, A_switch, cf, cb, c0, hf, hb, h0, sf, sb, Vaf, Vab, daf, dab)

        #-----------------------------------------------------------------
        # Collect/calculate Residuals
        #-----------------------------------------------------------------
        hjb_residual = vec(u) + A * vec(V) + V_dot + VEErrors - ρ_ss * vec(V)
        pc_residual  = -((r - 0) * inflation - (ceselast / priceadjust *
                                                (w / TFP - (ceselast-1) / ceselast) +
                                                inflation_dot - inflation_error))

        g_azdelta            = vec(g) .* vec(azdelta)
        g_intermediate       = spdiagm(0 => 1 ./ azdelta) * A' * g_azdelta
        g_residual           = g_dot - g_intermediate[1:end-1]

        mp_residual          = mp_dot - (-θ_MP * MP + σ_MP * mp_shock)

        realsav              = sum(vec(aa) .* vec(g) .* vec(azdelta))
        realsav_dot          = sum(vec(s)  .* vec(g) .* vec(azdelta))
        bondmarket_residual  = realsav_dot/realsav + govbcrule_fixnomB * inflation

        labmarket_residual   = sum(vec(zz) .* vec(h) .* vec(g) .* vec(azdelta)) - hours
        consumption_residual = sum(vec(c) .* vec(g) .* vec(azdelta)) - consumption
        output_residual      = TFP * hours - output
        assets_residual      = assets - realsav

        # Return equilibrium conditions
        return [hjb_residual; pc_residual; g_residual; mp_residual; bondmarket_residual; labmarket_residual;
                consumption_residual; output_residual; assets_residual]
    end

    derivs = ForwardDiff.jacobian(get_residuals, x)

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