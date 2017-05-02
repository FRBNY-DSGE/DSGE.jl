"""
```
measurement{T<:AbstractFloat}(m::Model1010{T}, TTT::Matrix{T}, RRR::Matrix{T},
                              CCC::Matrix{T}; shocks::Bool = true)
```

Assign measurement equation
```
X_t = ZZ*S_t + DD + u_t
```
where
```
u_t = eta_t + MM*eps_t
var(eta_t) = EE
var(u_t) = HH = EE + MM*QQ*MM'
cov(eps_t,u_t) = VV = QQ*MM'
```
"""
function measurement{T<:AbstractFloat}(m::Model1010{T},
                                       TTT::Matrix{T},
                                       RRR::Matrix{T},
                                       CCC::Vector{T};
                                       shocks::Bool = true)
    endo = m.endogenous_states
    exo  = m.exogenous_shocks
    obs  = m.observables

    # If shocks = true, then return measurement equation matrices with rows and columns for
    # anticipated policy shocks
    if shocks
        _n_observables = n_observables(m)
        _n_states = n_states_augmented(m)
        _n_shocks_exogenous = n_shocks_exogenous(m)
        endo_new = m.endogenous_states_augmented
    else
        _n_observables = n_observables(m) - n_anticipated_shocks(m)
        _n_states = n_states_augmented(m) - n_anticipated_shocks(m)
        _n_shocks_exogenous = n_shocks_exogenous(m) - n_anticipated_shocks(m)
        endo_new = Dict(
            [(key,m.endogenous_states_augmented[key] - n_anticipated_shocks(m)) for key in keys(m.endogenous_states_augmented)])
    end

    ZZ = zeros(_n_observables, _n_states)
    DD = zeros(_n_observables)
    MM = zeros(_n_observables, _n_shocks_exogenous)
    EE = zeros(_n_observables, _n_observables)
    QQ = zeros(_n_shocks_exogenous, _n_shocks_exogenous)

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
    ZZ[obs[:obs_wages], endo[:w_t]]      = 1.0
    ZZ[obs[:obs_wages], endo_new[:w_t1]] = -1.0
    ZZ[obs[:obs_wages], endo[:z_t]]      = 1.0
    DD[obs[:obs_wages]]                  = 100*(exp(m[:z_star])-1)

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

    ## 20 yrs forward transition matrix
    TTT20 = if subspec(m) == "ss16"
        eye(size(TTT,1))
    else
        (1/80)*((UniformScaling(1.) - TTT) \ (UniformScaling(1.) - TTT^80))
    end
    ## BAA or BBB Spread
    ZZ[obs[:obs_BBBspread], endo[:ERtil_k_t]]    = 1.0
    ZZ[obs[:obs_BBBspread], endo[:R_t]]          = -1.0
    # the previous 2 lines create a row of ZZ which is a vector of
    # zeros except for the coefficients on the relevant states. We
    # multiply this by TTT20 to adjust for the 20-year maturity of the
    # bond, and reassign the result to the BBB spread row of the ZZ
    # matrix, overwriting the original coefficients. Then we add a
    # coefficient for the measurement error.
    ZZ[obs[:obs_BBBspread], :]                   = ZZ[obs[:obs_BBBspread], :]' * TTT20
    ZZ[obs[:obs_BBBspread], endo_new[:e_BBB_t]]  = 1.0
    DD[obs[:obs_BBBspread]]                      = 100*log(m[:spr]*m[:lnb_safe]*m[:lnb_liq])

    ## AAA Spread
    ZZ[obs[:obs_AAAspread], endo[:b_liq_t]]       = -(1.0-m[:λ_AAA])*(m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))/(1 - m[:h]*exp(-m[:z_star]))
    ZZ[obs[:obs_AAAspread], endo[:ERtil_k_t]]     = m[:λ_AAA]
    ZZ[obs[:obs_AAAspread], endo[:R_t]]           = -m[:λ_AAA]
    # the previous 3 lines create a row of ZZ which is a vector of
    # zeros except for the coefficients on the relevant states. We
    # multiply this by TTT20 to adjust for the 20-year maturity of the
    # bond, and reassign the result to the AAA spread row of the ZZ
    # matrix, overwriting the original coefficients. Then we add a
    # coefficient for the measurement error.
    ZZ[obs[:obs_AAAspread], :]                    = ZZ[obs[:obs_AAAspread], :]' * TTT20
    ZZ[obs[:obs_AAAspread], endo_new[:e_AAA_t]]   = 1.0
    DD[obs[:obs_AAAspread]]                       = 100*log(m[:lnb_liq]) + 100*m[:λ_AAA]*log(m[:spr]*m[:lnb_safe])

    ## 10 yrs infl exp
    TTT10                          = (1/40)*((UniformScaling(1.) - TTT)\(UniformScaling(1.) - TTT^40))
    ZZ[obs[:obs_longinflation], :] = TTT10[endo[:π_t], :]
    DD[obs[:obs_longinflation]]    = 100*(m[:π_star]-1)

    ## Long Rate
    ZZ[obs[:obs_longrate], :]               = ZZ[obs[:obs_nominalrate], :]' * TTT10
    ZZ[obs[:obs_longrate], endo_new[:lr_t]] = 1.0
    DD[obs[:obs_longrate]]                  = m[:Rstarn]

    ## TFP
    ZZ[obs[:obs_tfp], endo[:z_t]]       = (1-m[:α])*m[:Iendoα] + 1*(1-m[:Iendoα])
    ZZ[obs[:obs_tfp], endo_new[:tfp_t]] = 1.0
    ZZ[obs[:obs_tfp], endo[:u_t]]       = m[:α]/( (1-m[:α])*(1-m[:Iendoα]) + 1*m[:Iendoα] )
    ZZ[obs[:obs_tfp], endo_new[:u_t1]]  = -(m[:α]/( (1-m[:α])*(1-m[:Iendoα]) + 1*m[:Iendoα]) )

    QQ[exo[:g_sh], exo[:g_sh]]                  = m[:σ_g]^2
    QQ[exo[:b_liqtil_sh], exo[:b_liqtil_sh]]    = m[:σ_b_liqtil]^2
    QQ[exo[:b_liqp_sh], exo[:b_liqp_sh]]        = m[:σ_b_liqp]^2
    QQ[exo[:b_safetil_sh], exo[:b_safetil_sh]]  = m[:σ_b_safetil]^2
    QQ[exo[:b_safep_sh], exo[:b_safep_sh]]      = m[:σ_b_safep]^2
    QQ[exo[:μ_sh], exo[:μ_sh]]                  = m[:σ_μ]^2
    QQ[exo[:z_sh], exo[:z_sh]]                  = m[:σ_z]^2
    QQ[exo[:λ_f_sh], exo[:λ_f_sh]]              = m[:σ_λ_f]^2
    QQ[exo[:λ_w_sh], exo[:λ_w_sh]]              = m[:σ_λ_w]^2
    QQ[exo[:rm_sh], exo[:rm_sh]]                = m[:σ_r_m]^2
    QQ[exo[:σ_ω_sh], exo[:σ_ω_sh]]              = m[:σ_σ_ω]^2
    QQ[exo[:μ_e_sh], exo[:μ_e_sh]]              = m[:σ_μ_e]^2
    QQ[exo[:γ_sh], exo[:γ_sh]]                  = m[:σ_γ]^2
    QQ[exo[:π_star_sh], exo[:π_star_sh]]        = m[:σ_π_star]^2
    QQ[exo[:lr_sh], exo[:lr_sh]]                = m[:σ_lr]^2
    QQ[exo[:zp_sh], exo[:zp_sh]]                = m[:σ_z_p]^2
    QQ[exo[:tfp_sh], exo[:tfp_sh]]              = m[:σ_tfp]^2
    QQ[exo[:gdpdef_sh], exo[:gdpdef_sh]]        = m[:σ_gdpdef]^2
    QQ[exo[:corepce_sh], exo[:corepce_sh]]      = m[:σ_corepce]^2
    QQ[exo[:gdp_sh], exo[:gdp_sh]]              = m[:σ_gdp]^2
    QQ[exo[:gdi_sh], exo[:gdi_sh]]              = m[:σ_gdi]^2
    QQ[exo[:AAA_sh], exo[:AAA_sh]]              = m[:σ_AAA]^2
    QQ[exo[:BBB_sh], exo[:BBB_sh]]              = m[:σ_BBB]^2

    # These lines set the standard deviations for the anticipated shocks. They
    # are here no longer calibrated to the std dev of contemporaneous shocks,
    # as we had in 904
    if shocks
        for i = 1:n_anticipated_shocks(m)
            ZZ[obs[Symbol("obs_nominalrate$i")], :]              = ZZ[obs[:obs_nominalrate], :]' * (TTT^i)
            DD[obs[Symbol("obs_nominalrate$i")]]                 = m[:Rstarn]
            QQ[exo[Symbol("rm_shl$i")], exo[Symbol("rm_shl$i")]] = m[Symbol("σ_r_m$i")]^2
        end
    end

    # Adjustment to DD because measurement equation assumes CCC is the zero vector
    if any(CCC != 0)
        DD += ZZ*((UniformScaling(1) - TTT)\CCC)
    end

    HH    = EE + MM*QQ*MM'
    VV    = QQ*MM'
    VVall = [[RRR*QQ*RRR' RRR*VV];
             [VV'*RRR'    HH]]

    return Measurement(ZZ, DD, QQ, EE, MM, VVall)
end
