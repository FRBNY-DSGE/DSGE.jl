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
                     CCC::Vector{T}; reg::Int = 1) where {T<:AbstractFloat}

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
            ModelConstructors.toggle_regime!(para, reg)
        end
    end

    no_integ_inds = inds_states_no_integ_series(m)
    if haskey(m.endogenous_states, :pgap_t)
        no_integ_inds = setdiff(no_integ_inds, [m.endogenous_states[:pgap_t]])
    end
    if haskey(m.endogenous_states, :ygap_t)
        no_integ_inds = setdiff(no_integ_inds, [m.endogenous_states[:ygap_t]])
    end
    if haskey(m.endogenous_states, :rw_t)
        no_integ_inds = setdiff(no_integ_inds, [m.endogenous_states[:rw_t]])
    end
    if haskey(m.endogenous_states, :Rref_t)
        no_integ_inds = setdiff(no_integ_inds, [m.endogenous_states[:Rref_t]])
    end
    if (get_setting(m, :add_laborproductivity_measurement) || get_setting(m, :add_nominalgdp_level) ||
        get_setting(m, :add_cumulative)) || (haskey(m.endogenous_states, :pgap_t)) || (haskey(m.endogenous_states, :ygap_t)) ||
        (haskey(m.endogenous_states, :rw_t)) || (haskey(m.endogenous_states, :Rref_t))
        # Remove integrated states (e.g. states w/unit roots)
        TTT = @view TTT[no_integ_inds, no_integ_inds]
        CCC = @view CCC[no_integ_inds]
    end

    ## GDP growth - Quarterly!
    ZZ[obs[:obs_gdp], endo[:y_t]]          = 1.0
    ZZ[obs[:obs_gdp], endo_new[:y_t1]]     = -1.0
    ZZ[obs[:obs_gdp], endo[:z_t]]          = 1.0
    ZZ[obs[:obs_gdp], endo_new[:e_gdp_t]]  = 1.0
    ZZ[obs[:obs_gdp], endo_new[:e_gdp_t1]] = -m[:me_level]
    DD[obs[:obs_gdp]]                      = 100*(exp(m[:z_star])-1)

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

    TTT10                                      = (1/40) * ((Matrix{Float64}(I, size(TTT, 1), size(TTT,1))
                                                            - TTT) \ (TTT - TTT^41))
    ZZ[obs[:obs_longinflation], no_integ_inds] = TTT10[endo[:π_t], :]
    DD[obs[:obs_longinflation]]                = 100*(m[:π_star]-1)

    ## Long Rate
    ZZ[obs[:obs_longrate], no_integ_inds]     = ZZ[6, no_integ_inds]' * TTT10
    ZZ[obs[:obs_longrate], endo_new[:e_lr_t]] = 1.0
    DD[obs[:obs_longrate]]                    = m[:Rstarn]

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

    ## Set up structural shocks covariance matrix
    QQ[exo[:g_sh], exo[:g_sh]]            = m[:σ_g]^2
    QQ[exo[:b_sh], exo[:b_sh]]            = m[:σ_b]^2
    QQ[exo[:μ_sh], exo[:μ_sh]]            = m[:σ_μ]^2
    QQ[exo[:ztil_sh], exo[:ztil_sh]]      = m[:σ_ztil]^2
    QQ[exo[:λ_f_sh], exo[:λ_f_sh]]        = m[:σ_λ_f]^2
    QQ[exo[:λ_w_sh], exo[:λ_w_sh]]        = m[:σ_λ_w]^2
    QQ[exo[:rm_sh], exo[:rm_sh]]          = m[:σ_r_m]^2
    QQ[exo[:σ_ω_sh], exo[:σ_ω_sh]]        = m[:σ_σ_ω]^2
    QQ[exo[:μ_e_sh], exo[:μ_e_sh]]        = m[:σ_μ_e]^2
    QQ[exo[:γ_sh], exo[:γ_sh]]            = m[:σ_γ]^2
    QQ[exo[:π_star_sh], exo[:π_star_sh]]  = m[:σ_π_star]^2
    QQ[exo[:lr_sh], exo[:lr_sh]]          = m[:σ_lr]^2
    QQ[exo[:zp_sh], exo[:zp_sh]]          = m[:σ_z_p]^2
    QQ[exo[:tfp_sh], exo[:tfp_sh]]        = m[:σ_tfp]^2
    QQ[exo[:gdpdef_sh], exo[:gdpdef_sh]]  = m[:σ_gdpdef]^2
    QQ[exo[:corepce_sh], exo[:corepce_sh]]= m[:σ_corepce]^2
    QQ[exo[:gdp_sh], exo[:gdp_sh]]        = m[:σ_gdp]^2
    QQ[exo[:gdi_sh], exo[:gdi_sh]]        = m[:σ_gdi]^2

    if subspec(m) in ["ss59", "ss60", "ss61"]
        QQ[exo[:ziid_sh], exo[:ziid_sh]] = m[:σ_ziid]^2
        QQ[exo[:biidc_sh], exo[:biidc_sh]] = m[:σ_biidc]^2
        QQ[exo[:φ_sh], exo[:φ_sh]] = m[:σ_φ]^2
    end

    # Automated addition of anticipated shocks to QQ
    for (k, v) in get_setting(m, :antshocks)
        for i = 1:v
            QQ[exo[Symbol(k, "_shl$i")], exo[Symbol(k, "_shl$i")]] = m[Symbol("σ_", k, "$i")]^2
        end
    end

    ## Anticipated observables

    # Anticipated monetary policy shocks
    for i = 1:n_mon_anticipated_shocks(m)
        ZZ[obs[Symbol("obs_nominalrate$i")], no_integ_inds] = ZZ[obs[:obs_nominalrate], no_integ_inds]' * (TTT^i)
        DD[obs[Symbol("obs_nominalrate$i")]]    = m[:Rstarn]
        if subspec(m) == "ss11"
            QQ[exo[Symbol("rm_shl$i")], exo[Symbol("rm_shl$i")]] = m[:σ_r_m]^2 / n_mon_anticipated_shocks(m)
        else
            QQ[exo[Symbol("rm_shl$i")], exo[Symbol("rm_shl$i")]] = m[Symbol("σ_r_m$i")]^2
        end
    end

    # Anticipated GDP growth
    if haskey(get_settings(m), :add_anticipated_obs_gdp)
        if get_setting(m, :add_anticipated_obs_gdp)
            for i = 1:get_setting(m, :n_anticipated_obs_gdp)
                ZZ_obs_gdp = ZZ[obs[:obs_gdp], :]
                ZZ_obs_gdp[endo_new[:e_gdp_t]]  = 0. # Ignore measurement error for anticipated GDP growth
                ZZ_obs_gdp[endo_new[:e_gdp_t1]] = 0.
                ZZ[obs[Symbol("obs_gdp$i")], no_integ_inds] = ZZ_obs_gdp[no_integ_inds]' * (TTT^i)
                DD[obs[Symbol("obs_gdp$i")]]                = 100. * (exp(m[:z_star]) - 1.)
            end
        end
    end

    # Adjustment to DD because measurement equation assumes CCC is the zero vector
    if any(CCC .!= 0.)
        # Are we using the zero rate alternative policy or not?
        is_zero_rate_rule = (haskey(get_settings(m), :replace_eqcond_func_dict) ?
                             (!isempty(get_setting(m, :replace_eqcond_func_dict)) && get_setting(m, :replace_eqcond)) : false) ||
                             get_setting(m, :alternative_policy).eqcond == zero_rate_eqcond

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
            DD += ZZ[:, no_integ_inds]*((UniformScaling(1) - TTT) \ CCC)
        end
    end

    for para in m.parameters
        if !isempty(para.regimes)
            ModelConstructors.toggle_regime!(para, 1)
        end
    end

    return Measurement(ZZ, DD, QQ, EE)
end
