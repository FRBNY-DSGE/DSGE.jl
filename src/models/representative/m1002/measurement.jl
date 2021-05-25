"""
```
measurement(m::Model1002{T}, TTT::Matrix{T}, RRR::Matrix{T},
    CCC::Vector{T}; reg::Int = 1) where {T<:AbstractFloat}
```

Assign measurement equation

```
y_t = ZZ*s_t + DD + u_t
```

where

```
Var(ϵ_t) = QQ
Var(u_t) = EE
Cov(ϵ_t, u_t) = 0
```
"""
function measurement(m::Model1002{T},
                     TTT::Matrix{T},
                     RRR::Matrix{T},
                     CCC::Vector{T}; reg::Int = 1,
                     TTTs::Vector{<: AbstractMatrix{T}} = Matrix{T}[],
                     CCCs::Vector{<: AbstractVector{T}} = Vector{T}[],
                     information_set::UnitRange = reg:reg,
                     memo::Union{ForwardMultipleExpectationsMemo, Nothing} = nothing) where {T <: AbstractFloat}

    endo     = m.endogenous_states
    endo_new = m.endogenous_states_augmented
    exo      = m.exogenous_shocks
    obs      = m.observables

    _n_observables = n_observables(m)
    _n_states = n_states_augmented(m)
    _n_shocks_exogenous = n_shocks_exogenous(m)

    ZZ = zeros(_n_observables, _n_states)
    DD = zeros(_n_observables)
    EE = zeros(_n_observables, _n_observables)
    QQ = zeros(_n_shocks_exogenous, _n_shocks_exogenous)

    for para in m.parameters
        if !isempty(para.regimes)
            if haskey(get_settings(m), :model2para_regime) && haskey(get_setting(m, :model2para_regime), para.key)
                ModelConstructors.toggle_regime!(para, reg, get_setting(m, :model2para_regime)[para.key])
            else
                ModelConstructors.toggle_regime!(para, reg)
            end
        end
    end

    # Set up for calculating k-periods ahead expectations and expected sums
    permanent_t = length(information_set[findfirst(x -> x == reg, information_set):end]) - 1 + reg
    if information_set[1] == information_set[end] # this condition is probably wrong
        # In this case, we do not need to pass the TTTs, CCCs in, so we redefine them as empty.
        # This step is also necessary to ensure the memo is properly used.
        TTTs = Matrix{T}[]
        CCCs = Vector{T}[]

        # for now, assume that in this case, we do not want to use the memo.
        # The complicating problem is that, in the historical regimes, we don't want to use the memo
        # (at least, not the one that is currently being automatically constructed).
        # But in the forecast regimes, we do want to use the memo if we ever reach a case where
        # the information set becomes myopic.
        memo = nothing
        # TODO: maybe instead of emptying TTTs, CCCs, we add
        # a step to recompute the memo since we will still likely be recalculating
        # products/powers of TTT multiple times that could be pre-computed.
    end
    use_fwd_exp_sum = haskey(get_settings(m), :use_forward_expected_sum_memo) && get_setting(m, :use_forward_expected_sum_memo)
    use_fwd_exp     = haskey(get_settings(m), :use_forward_expectations_memo) && get_setting(m, :use_forward_expectations_memo)

    flex_ait_2020Q3 = haskey(get_settings(m), :flexible_ait_policy_change) &&
        get_setting(m, :flexible_ait_policy_change)

    # With flexible_ait_policy_change: there is a regime-break in 2020:Q3, so before 2020:Q2,
    # the final regime should be considered to be the regime corresponding to 2020:Q2, namely for regimes
    # before 2020:Q3, the measurement equation should reflect the belief that the "permanent" regime is the
    # regime before 2020:Q3 starts.
    if flex_ait_2020Q3
        if get_setting(m, :regime_dates)[4] != get_setting(m, :flexible_ait_policy_change_date)
            reg_2020Q3 = 0
            for i in 1:get_setting(m, :n_regimes)
                if get_setting(m, :regime_dates)[i] == get_setting(m, :flexible_ait_policy_change_date)
                    reg_2020Q2 = i
                    break
                end
            end
            @assert reg_2020Q2 != 0 "The setting :regime_dates does not contain the date 2020:Q3, which is required if the setting :flexible_ait_policy_change is true."
        else
            reg_2020Q3 = 4
        end

        permanent_t = reg < reg_2020Q3 ? reg_2020Q3 - 1 : permanent_t
    end

    # Remove integrated states (e.g. states w/unit roots)
    no_integ_inds = inds_states_no_integ_series(m)
    if ((haskey(get_settings(m), :add_altpolicy_pgap) && haskey(get_settings(m), :pgap_type)) ?
        (get_setting(m, :add_altpolicy_pgap) && get_setting(m, :pgap_type) == :ngdp) : false) ||
        (haskey(get_settings(m), :ait_Thalf) ? (exp(log(0.5) / get_setting(m, :ait_Thalf)) ≈ 0.) : false)
        no_integ_inds = setdiff(no_integ_inds, [m.endogenous_states[:pgap_t]])
    end
    if (haskey(get_settings(m), :gdp_Thalf) ? (exp(log(0.5) / get_setting(m, :gdp_Thalf)) ≈ 0.) : false)
        no_integ_inds = setdiff(no_integ_inds, [m.endogenous_states[:ygap_t]])
    end
    if haskey(get_settings(m), :ρ_rw) ? (get_setting(m, :ρ_rw) ≈ 1.) : false
        no_integ_inds = setdiff(no_integ_inds, [m.endogenous_states[:rw_t]])
    end
    if haskey(get_settings(m), :rw_ρ_smooth) ? (get_setting(m, :rw_ρ_smooth) ≈ 1.) : false
        no_integ_inds = setdiff(no_integ_inds, [m.endogenous_states[:Rref_t]])
    end
    integ_series = length(no_integ_inds) != n_states_augmented(m) # Are the series integrated?

    ## GDP growth - Quarterly!
    ZZ[obs[:obs_gdp], endo[:y_t]]          = 1.0
    ZZ[obs[:obs_gdp], endo_new[:y_t1]]     = -1.0
    ZZ[obs[:obs_gdp], endo[:z_t]]          = 1.0
    ZZ[obs[:obs_gdp], endo_new[:e_gdp_t]]  = 1.0
    ZZ[obs[:obs_gdp], endo_new[:e_gdp_t1]] = -m[:me_level]
    DD[obs[:obs_gdp]]                      = 100*(exp(m[:z_star])-1)

    if haskey(get_settings(m), :add_iid_cond_obs_gdp_meas_err) ?
        get_setting(m, :add_iid_cond_obs_gdp_meas_err) : false
        ZZ[obs[:obs_gdp], endo_new[:e_condgdp_t]] = 1.
    end

    ## GDI growth- Quarterly!
    ZZ[obs[:obs_gdi], endo[:y_t]]          = m[:γ_gdi]
    ZZ[obs[:obs_gdi], endo_new[:y_t1]]     = -m[:γ_gdi]
    ZZ[obs[:obs_gdi], endo[:z_t]]          = m[:γ_gdi]
    ZZ[obs[:obs_gdi], endo_new[:e_gdi_t]]  = 1.0
    ZZ[obs[:obs_gdi], endo_new[:e_gdi_t1]] = -m[:me_level]
    DD[obs[:obs_gdi]]                      = 100*(exp(m[:z_star])-1) + m[:δ_gdi]

    ## Hours growth
    ZZ[obs[:obs_hours], endo[:L_t]] = 1.0
    DD[obs[:obs_hours]]             = m[:Lmean]

    ## Labor Share/real wage growth
    if subspec(m) in ["ss16", "ss17"]
        # log(labor_share) = log(wage) + log(hours) - log(GDP)
        ZZ[obs[:obs_laborshare], endo[:w_t]] = 1.
        ZZ[obs[:obs_laborshare], endo[:L_t]] = 1.
        DD[obs[:obs_laborshare]] = 100. * log(m[:wstar] * m[:Lstar] / m[:ystar])
        ZZ[obs[:obs_laborshare], endo[:y_t]] = -1.
    else
        ZZ[obs[:obs_wages], endo[:w_t]]      = 1.0
        ZZ[obs[:obs_wages], endo_new[:w_t1]] = -1.0
        ZZ[obs[:obs_wages], endo[:z_t]]      = 1.0
        DD[obs[:obs_wages]]                  = 100*(exp(m[:z_star])-1)
    end

    ## Inflation (GDP Deflator)
    ZZ[obs[:obs_gdpdeflator], endo[:π_t]]            = m[:Γ_gdpdef]
    ZZ[obs[:obs_gdpdeflator], endo_new[:e_gdpdef_t]] = 1.0
    DD[obs[:obs_gdpdeflator]]                        = 100*(m[:π_star]-1) + m[:δ_gdpdef]

    ## Inflation (Core PCE)
    ZZ[obs[:obs_corepce], endo[:π_t]]             = 1.0
    ZZ[obs[:obs_corepce], endo_new[:e_corepce_t]] = 1.0
    DD[obs[:obs_corepce]]                         = 100*(m[:π_star]-1)

    if haskey(get_settings(m), :add_iid_cond_obs_corepce_meas_err) ?
        get_setting(m, :add_iid_cond_obs_corepce_meas_err) : false
        ZZ[obs[:obs_corepce], endo_new[:e_condcorepce_t]] = 1.
    end

    ## Nominal interest rate
    ZZ[obs[:obs_nominalrate], endo[:R_t]] = 1.0
    DD[obs[:obs_nominalrate]]             = m[:Rstarn]

    ## Consumption Growth
    ZZ[obs[:obs_consumption], endo[:c_t]]      = 1.0
    ZZ[obs[:obs_consumption], endo_new[:c_t1]] = -1.0
    ZZ[obs[:obs_consumption], endo[:z_t]]      = 1.0
    DD[obs[:obs_consumption]]                  = 100*(exp(m[:z_star])-1)

    ## Investment Growth
    ZZ[obs[:obs_investment], endo[:i_t]]      = 1.0
    ZZ[obs[:obs_investment], endo_new[:i_t1]] = -1.0
    ZZ[obs[:obs_investment], endo[:z_t]]      = 1.0
    DD[obs[:obs_investment]]                  = 100*(exp(m[:z_star])-1)

    ## Spreads
    ZZ[obs[:obs_spread], endo[:ERktil_t]] = 1.0
    ZZ[obs[:obs_spread], endo[:R_t]]       = -1.0
    DD[obs[:obs_spread]]                   = 100*log(m[:spr])

    ## 10 yrs infl exp
    TTT10, CCC10 = k_periods_ahead_expected_sums(TTT, CCC, TTTs, CCCs, reg, 40, permanent_t;
                                                 integ_series = integ_series,
                                                 memo = use_fwd_exp_sum ? memo : nothing)
    TTT10        = TTT10 ./ 40. # divide by 40 to average across 10 years
    CCC10        = CCC10 ./ 40.
    ZZ[obs[:obs_longinflation], :] = view(TTT10, endo[:π_t], :)
    # ZZ[obs[:obs_longinflation], :] = TTT10[endo[:π_t], :]
    DD[obs[:obs_longinflation]]    = 100*(m[:π_star]-1) + CCC10[endo[:π_t]]

    ## Long Rate
    # ZZ_long_rate = ZZ[obs[:obs_nominalrate], :]'
    # ZZ[obs[:obs_longrate], :]                 = ZZ_long_rate * TTT10 # TODO: this is slightly inefficient, do what we do with long inflation
    ZZ[obs[:obs_longrate], :]                 =  view(TTT10, endo[:R_t], :) # TODO: this is slightly inefficient, do what we do with long inflation
    ZZ[obs[:obs_longrate], endo_new[:e_lr_t]] = 1.0
    DD[obs[:obs_longrate]]                    = m[:Rstarn] + CCC10[endo[:R_t]]
    # DD[obs[:obs_longrate]]                    = m[:Rstarn] + ZZ_long_rate * CCC10

    ## TFP
    ZZ[obs[:obs_tfp], endo[:z_t]] = (1-m[:α])*m[:Iendoα] + 1*(1-m[:Iendoα])
    if subspec(m) in ["ss14", "ss15", "ss16", "ss18", "ss19"]
        ZZ[obs[:obs_tfp], endo_new[:e_tfp_t]]  = 1.0
        ZZ[obs[:obs_tfp], endo_new[:e_tfp_t1]] = -m[:me_level]
    else
        ZZ[obs[:obs_tfp], endo_new[:e_tfp_t]] = 1.0
    end
    if !(subspec(m) in ["ss15", "ss16"])
        ZZ[obs[:obs_tfp], endo[:u_t]]       = m[:α]/( (1-m[:α])*(1-m[:Iendoα]) + 1*m[:Iendoα] )
        ZZ[obs[:obs_tfp], endo_new[:u_t1]]  = -(m[:α]/( (1-m[:α])*(1-m[:Iendoα]) + 1*m[:Iendoα]) )
    end

    # ygap and pgap for Flexible AIT rule
    if (haskey(get_settings(m), :add_initialize_pgap_ygap_pseudoobs) ? get_setting(m, :add_initialize_pgap_ygap_pseudoobs) : false)
        ZZ[obs[:obs_pgap], endo[:pgap_t]] = 1.
        ZZ[obs[:obs_ygap], endo[:ygap_t]] = 1.
    end

    ## Set up structural shocks covariance matrix
    QQ[exo[:g_sh], exo[:g_sh]]             = m[:σ_g]^2
    QQ[exo[:b_sh], exo[:b_sh]]             = m[:σ_b]^2
    QQ[exo[:μ_sh], exo[:μ_sh]]             = m[:σ_μ]^2
    QQ[exo[:ztil_sh], exo[:ztil_sh]]       = m[:σ_ztil]^2
    QQ[exo[:λ_f_sh], exo[:λ_f_sh]]         = m[:σ_λ_f]^2
    QQ[exo[:λ_w_sh], exo[:λ_w_sh]]         = m[:σ_λ_w]^2
    QQ[exo[:rm_sh], exo[:rm_sh]]           = m[:σ_r_m]^2
    QQ[exo[:σ_ω_sh], exo[:σ_ω_sh]]         = m[:σ_σ_ω]^2
    QQ[exo[:μ_e_sh], exo[:μ_e_sh]]         = m[:σ_μ_e]^2
    QQ[exo[:γ_sh], exo[:γ_sh]]             = m[:σ_γ]^2
    QQ[exo[:π_star_sh], exo[:π_star_sh]]   = m[:σ_π_star]^2
    QQ[exo[:lr_sh], exo[:lr_sh]]           = m[:σ_lr]^2
    QQ[exo[:zp_sh], exo[:zp_sh]]           = m[:σ_z_p]^2
    QQ[exo[:tfp_sh], exo[:tfp_sh]]         = m[:σ_tfp]^2
    QQ[exo[:gdpdef_sh], exo[:gdpdef_sh]]   = m[:σ_gdpdef]^2
    QQ[exo[:corepce_sh], exo[:corepce_sh]] = m[:σ_corepce]^2
    QQ[exo[:gdp_sh], exo[:gdp_sh]]         = m[:σ_gdp]^2
    QQ[exo[:gdi_sh], exo[:gdi_sh]]         = m[:σ_gdi]^2

    if subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63", "ss64", "ss65", "ss66", "ss67", "ss68", "ss69", "ss70", "ss71", "ss72", "ss73", "ss74", "ss75", "ss76", "ss77", "ss78", "ss79", "ss80", "ss81", "ss82", "ss83"]
        QQ[exo[:ziid_sh], exo[:ziid_sh]]   = m[:σ_ziid]^2
        QQ[exo[:biidc_sh], exo[:biidc_sh]] = m[:σ_biidc]^2
        QQ[exo[:φ_sh], exo[:φ_sh]]         = m[:σ_φ]^2
    end

    if subspec(m) in ["ss67", "ss68", "ss69", "ss70", "ss71", "ss72", "ss73", "ss74", "ss75", "ss76", "ss77", "ss78", "ss80", "ss82", "ss83"]
        QQ[exo[:g_covid_sh], exo[:g_covid_sh]]   = m[:σ_g_covid]^2
        QQ[exo[:μ_covid_sh], exo[:μ_covid_sh]]   = m[:σ_μ_covid]^2
        QQ[exo[:λ_f_covid_sh], exo[:λ_f_covid_sh]]   = m[:σ_λ_f_covid]^2
        QQ[exo[:σ_ω_covid_sh], exo[:σ_ω_covid_sh]]   = m[:σ_σ_ω_covid]^2
        QQ[exo[:lr_covid_sh], exo[:lr_covid_sh]]   = m[:σ_lr_covid]^2
        QQ[exo[:tfp_covid_sh], exo[:tfp_covid_sh]]   = m[:σ_tfp_covid]^2
        QQ[exo[:gdp_covid_sh], exo[:gdp_covid_sh]]   = m[:σ_gdp_covid]^2
        QQ[exo[:gdi_covid_sh], exo[:gdi_covid_sh]]   = m[:σ_gdi_covid]^2
    end
    if subspec(m) in ["ss69", "ss70", "ss73", "ss74", "ss77", "ss78"]
        QQ[exo[:zp_covid_sh], exo[:zp_covid_sh]]   = m[:σ_z_p_covid]^2
    end

    if subspec(m) in ["ss62"] # TODO: add ss63
        ## Set up structural shocks covariance matrix
        QQ[exo[:g_sh], exo[:g_sh]]             *= m[:damp_standard_shocks]^2
        QQ[exo[:b_sh], exo[:b_sh]]             *= m[:damp_standard_shocks]^2
        QQ[exo[:μ_sh], exo[:μ_sh]]             *= m[:damp_standard_shocks]^2
        QQ[exo[:ztil_sh], exo[:ztil_sh]]       *= m[:damp_standard_shocks]^2
        QQ[exo[:λ_f_sh], exo[:λ_f_sh]]         *= m[:damp_standard_shocks]^2
        QQ[exo[:λ_w_sh], exo[:λ_w_sh]]         *= m[:damp_standard_shocks]^2
        QQ[exo[:rm_sh], exo[:rm_sh]]           *= m[:amplify_σ_r_m]^2
        QQ[exo[:σ_ω_sh], exo[:σ_ω_sh]]         *= m[:damp_standard_shocks]^2
        QQ[exo[:μ_e_sh], exo[:μ_e_sh]]         *= m[:damp_standard_shocks]^2
        QQ[exo[:γ_sh], exo[:γ_sh]]             *= m[:damp_standard_shocks]^2
        QQ[exo[:π_star_sh], exo[:π_star_sh]]   *= m[:damp_standard_shocks]^2
        QQ[exo[:lr_sh], exo[:lr_sh]]           *= m[:damp_standard_shocks]^2
        # QQ[exo[:zp_sh], exo[:zp_sh]]           = m[:σ_z_p]^2 # σ_z_p is handled on its own via regime-switching
        QQ[exo[:tfp_sh], exo[:tfp_sh]]         *= m[:damp_standard_shocks]^2
        QQ[exo[:gdpdef_sh], exo[:gdpdef_sh]]   *= m[:amplify_inflation_me]^2
        QQ[exo[:corepce_sh], exo[:corepce_sh]] *= m[:amplify_inflation_me]^2
        QQ[exo[:gdp_sh], exo[:gdp_sh]]         *= m[:damp_standard_shocks]^2
        QQ[exo[:gdi_sh], exo[:gdi_sh]]         *= m[:damp_standard_shocks]^2
    end

    # ygap and pgap shocks
    if (haskey(get_settings(m), :add_initialize_pgap_ygap_pseudoobs) ? get_setting(m, :add_initialize_pgap_ygap_pseudoobs) : false)
        QQ[exo[:pgap_sh], exo[:pgap_sh]] = m[:σ_pgap]^2
        QQ[exo[:ygap_sh], exo[:ygap_sh]] = m[:σ_ygap]^2
    end

    # measurement errors for "conditional" observations of GDP
    if haskey(get_settings(m), :add_iid_cond_obs_gdp_meas_err) ?
        get_setting(m, :add_iid_cond_obs_gdp_meas_err) : false
        QQ[exo[:condgdp_sh], exo[:condgdp_sh]] = m[:σ_condgdp] ^ 2
    end

    if haskey(get_settings(m), :add_iid_anticipated_obs_gdp_meas_err) ?
        get_setting(m, :add_iid_anticipated_obs_gdp_meas_err) : false
        QQ[exo[:gdpexp_sh], exo[:gdpexp_sh]] = m[:σ_gdpexp] ^ 2
    end

    # measurement errors for "conditional" observations of Core PCE
    if haskey(get_settings(m), :add_iid_cond_obs_corepce_meas_err) ?
        get_setting(m, :add_iid_cond_obs_corepce_meas_err) : false
        QQ[exo[:condcorepce_sh], exo[:condcorepce_sh]] = m[:σ_condcorepce] ^ 2
    end

    # Automated addition of anticipated shocks to QQ
    for (k, v) in get_setting(m, :antshocks)
        for i = 1:v
            QQ[exo[Symbol(k, "_shl$i")], exo[Symbol(k, "_shl$i")]] = m[Symbol("σ_", k, "$i")]^2
        end
    end

    ## Anticipated observables
    use_current_regime = haskey(get_settings(m), :measurement_use_current_regime_matrices) ?
        get_setting(m, :measurement_use_current_regime_matrices) : true

    # Anticipated monetary policy shocks
    # ZZ_obs_nomrate = ZZ[obs[:obs_nominalrate], :]'
    for i = 1:n_mon_anticipated_shocks(m)
        TTT_accum, CCC_accum = k_periods_ahead_expectations(TTT, CCC, TTTs, CCCs, reg, i, permanent_t;
                                                            integ_series = integ_series,
                                                            memo = (isnothing(memo) || !use_fwd_exp) ? nothing :
                                                            ForwardExpectationsMemo(memo.time_varying_memo[min(reg + i, permanent_t, 17)],
                                                                                    memo.permanent_memo))

        # ZZ[obs[Symbol("obs_nominalrate$i")], :] = ZZ_obs_nomrate * TTT_accum
        ZZ[obs[Symbol("obs_nominalrate$i")], :] = view(TTT_accum, endo[:R_t], :)
        # DD[obs[Symbol("obs_nominalrate$i")]]    = m[:Rstarn] + ZZ_obs_nomrate * CCC_accum
        DD[obs[Symbol("obs_nominalrate$i")]]    = m[:Rstarn] + CCC_accum[endo[:R_t]]
        if subspec(m) == "ss11"
            QQ[exo[Symbol("rm_shl$i")], exo[Symbol("rm_shl$i")]] = m[:σ_r_m]^2 / n_mon_anticipated_shocks(m)
        else
            QQ[exo[Symbol("rm_shl$i")], exo[Symbol("rm_shl$i")]] = m[Symbol("σ_r_m$i")]^2
        end
    end

    # Anticipated GDP growth
    if haskey(get_settings(m), :add_anticipated_obs_gdp)
        if get_setting(m, :add_anticipated_obs_gdp)
            ZZ_obs_gdp = ZZ[obs[:obs_gdp], :]'
            meas_err = haskey(get_settings(m), :meas_err_anticipated_obs_gdp) ?
                get_setting(m, :meas_err_anticipated_obs_gdp) : 0. # Ignore measurement error for anticipated GDP growth
            ZZ_obs_gdp[endo_new[:e_gdp_t]]  = meas_err
            ZZ_obs_gdp[endo_new[:e_gdp_t1]] = -meas_err * m[:me_level]

            for i = 1:get_setting(m, :n_anticipated_obs_gdp)
                TTT_accum, CCC_accum = k_periods_ahead_expectations(TTT, CCC, TTTs, CCCs, reg, i, permanent_t;
                                                                    integ_series = integ_series,
                                                                    memo = (isnothing(memo) || !use_fwd_exp) ? nothing :
                                                                    ForwardExpectationsMemo(memo.time_varying_memo[min(reg + i, permanent_t, 17)],
                                                                                            memo.permanent_memo))
                ZZ[obs[Symbol("obs_gdp$i")], :] = ZZ_obs_gdp * TTT_accum
                DD[obs[Symbol("obs_gdp$i")]]    = 100. * (exp(m[:z_star]) - 1.) + ZZ_obs_gdp * CCC_accum
                if haskey(get_settings(m), :add_iid_anticipated_obs_gdp_meas_err) ?
                    get_setting(m, :add_iid_anticipated_obs_gdp_meas_err) : false
                    ZZ[obs[Symbol("obs_gdp$i")], endo_new[:e_gdpexp_t]] = 1.
                end
            end
        end
    end

    # Adjustment to DD because measurement equation assumes CCC is the zero vector
    if any(CCC .!= 0.)
        # Are we using the zero rate alternative policy or not?
        is_zero_rate_rule = (haskey(get_settings(m), :regime_eqcond_info) ?
                             (!isempty(get_setting(m, :regime_eqcond_info)) && get_setting(m, :replace_eqcond)) : false) ||
        alternative_policy(m).eqcond == zero_rate_eqcond

        # If we are using the zero rate rule, then we don't adjust DD. If we run the code block below,
        # then the zero rate rule will cause CCC to be nonzero in multiple places, leading to the
        # wrong steady state measurement vector DD.
        #
        # The current coding of is_zero_rate_rule assumes that, when using a zero_rate rule (perhaps with gensys2 and temp alt policies)
        # any regimes in which the zero rate rule does not apply also will not have to update DD to reflect nonzero CCC.
        # This is typically the case for m1002.
        #
        # Currently it looks like we modify CCC in augment_states to track expected inflation,
        # but it doesn't look like we use this variable anywhere anyway. So ignoring this code block is fine on that front.
        # However, all the observables and pseudo-observables were not thoroughly checked, so it is possible
        # that the measurement equation is flawed in some way given how this has been coded.
        if !is_zero_rate_rule
            DD += ZZ * ((UniformScaling(1) - TTT) \ CCC)
        end
    end

    for para in m.parameters
        if !isempty(para.regimes)
            ModelConstructors.toggle_regime!(para, 1)
        end
    end

    return Measurement(ZZ, DD, QQ, EE)
end
