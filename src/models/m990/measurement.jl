# Assign measurement equation : X_t = ZZ*S_t + DD + u_t
# where u_t = eta_t+MM* eps_t with var(eta_t) = EE
# where var(u_t) = HH = EE+MM QQ MM', cov(eps_t,u_t) = VV = QQ*MM'

function measurement(m::Model990, TTT::Matrix, RRR::Matrix, CCC::Matrix; shocks::Bool = true)
    endo, exo, obs  = m.endogenous_states, m.exogenous_shocks, m.observables

    # If shocks = true, then return measurement equation matrices with rows and columns for anticipated policy shocks
    if shocks
        _observables = observables(m)
        _states = augmented_states(m)
        _exogenous_shocks = exogenous_shocks(m)
        endo_addl = m.endogenous_states_postgensys
    else
        _observables = observables(m) - anticipated_shocks(m)
        _states = augmented_states(m) - anticipated_shocks(m)
        _exogenous_shocks = exogenous_shocks(m) - anticipated_shocks(m)
        endo_addl = Dict(
            [(key,m.endogenous_states_postgensys[key] - anticipated_shocks(m)) for key in keys(m.endogenous_states_postgensys)])
    end

    ZZ = zeros(_observables, _states)
    DD = zeros(_observables, 1)
    MM = zeros(_observables, _exogenous_shocks)
    EE = zeros(_observables, _observables)
    QQ = zeros(_exogenous_shocks, _exogenous_shocks)

    ## Output growth - Quarterly!
    ZZ[obs[:g_y], endo[:y_t]]       = 1.0
    ZZ[obs[:g_y], endo_addl[:y_t1]] = -1.0
    ZZ[obs[:g_y], endo[:z_t]]       = 1.0
    DD[obs[:g_y]]                   = 100*(exp(m[:zstar])-1)

    ## Hoursg
    ZZ[obs[:hoursg], endo[:L_t]] = 1.0
    DD[obs[:hoursg]]             = m[:Lmean].scaledvalue

    ## Labor Share/real wage growth
    ZZ[obs[:g_w], endo[:w_t]]       = 1.0
    ZZ[obs[:g_w], endo_addl[:w_t1]] = -1.0
    ZZ[obs[:g_w], endo[:z_t]]       = 1.0
    DD[obs[:g_w]]                   = 100*(exp(m[:zstar])-1)

    ## Inflation (GDP Deflator)
    ZZ[obs[:pi_gdpdef], endo[:pi_t]]          = m[:gamm_gdpdef]
    ZZ[obs[:pi_gdpdef], endo_addl[:e_gdpdef]] = 1.0
    DD[obs[:pi_gdpdef]]                       = 100*(m[:pistar]-1) + m[:del_gdpdef]

    ## Inflation (Core PCE)
    ZZ[obs[:pi_pce], endo[:pi_t]]       = 1.0
    ZZ[obs[:pi_pce], endo_addl[:e_pce]] = 1.0
    DD[obs[:pi_pce]]                    = 100*(m[:pistar]-1)

    ## Nominal interest rate
    ZZ[obs[:R_n], endo[:R_t]] = 1.0
    DD[obs[:R_n]]              = m[:Rstarn]

    ## Consumption Growth
    ZZ[obs[:g_c], endo[:c_t]]       = 1.0
    ZZ[obs[:g_c], endo_addl[:c_t1]] = -1.0
    ZZ[obs[:g_c], endo[:z_t]]       = 1.0
    DD[obs[:g_c]]                   = 100*(exp(m[:zstar])-1)

    ## Investment Growth
    ZZ[obs[:g_i], endo[:i_t]]       = 1.0
    ZZ[obs[:g_i], endo_addl[:i_t1]] = -1.0
    ZZ[obs[:g_i], endo[:z_t]]       = 1.0
    DD[obs[:g_i]]                    = 100*(exp(m[:zstar])-1)

    ## Spreads
    ZZ[obs[:sprd], endo[:E_Rktil]] = 1.0
    ZZ[obs[:sprd], endo[:R_t]]     = -1.0
    DD[obs[:sprd]]                 = 100*log(m[:sprd])

    ## 10 yrs infl exp
    TTT10                = (1/40)*((eye(size(TTT, 1)) - TTT)\(TTT - TTT^41))
    ZZ[obs[:pi_long], :] =  TTT10[endo[:pi_t], :]
    DD[obs[:pi_long]]    = 100*(m[:pistar]-1)

    ## Long Rate
    ZZ[obs[:R_long], :]                = ZZ[6, :]*TTT10
    ZZ[obs[:R_long], endo_addl[:lr_t]] = 1.0
    DD[obs[:R_long]]                   = m[:Rstarn]

    ## TFP
    ZZ[obs[:tfp], endo[:z_t]]           = (1-m[:alp])*m[:modelalp_ind] + 1*(1-m[:modelalp_ind])
    ZZ[obs[:tfp], endo_addl[:tfp_t]]    = 1.0
    ZZ[obs[:tfp], endo[:u_t]]           = m[:alp]/( (1-m[:alp])*(1-m[:modelalp_ind]) + 1*m[:modelalp_ind] )
    ZZ[obs[:tfp], endo_addl[:u_t1]]     = -(m[:alp]/( (1-m[:alp])*(1-m[:modelalp_ind]) + 1*m[:modelalp_ind]) )

    QQ[exo[:g_sh], exo[:g_sh]]           = m[:σ_g]^2
    QQ[exo[:b_sh], exo[:b_sh]]           = m[:σ_b]^2
    QQ[exo[:mu_sh], exo[:mu_sh]]         = m[:σ_mu]^2
    QQ[exo[:z_sh], exo[:z_sh]]           = m[:σ_z]^2
    QQ[exo[:laf_sh], exo[:laf_sh]]       = m[:σ_laf]^2
    QQ[exo[:law_sh], exo[:law_sh]]       = m[:σ_law]^2
    QQ[exo[:rm_sh], exo[:rm_sh]]         = m[:σ_rm]^2
    QQ[exo[:sigw_sh], exo[:sigw_sh]]     = m[:σ_sigw]^2
    QQ[exo[:mue_sh], exo[:mue_sh]]       = m[:σ_mue]^2
    QQ[exo[:gamm_sh], exo[:gamm_sh]]     = m[:σ_gamm]^2
    QQ[exo[:pist_sh], exo[:pist_sh]]     = m[:σ_pist]^2
    QQ[exo[:lr_sh], exo[:lr_sh]]         = m[:σ_lr]^2
    QQ[exo[:zp_sh], exo[:zp_sh]]         = m[:σ_zp]^2
    QQ[exo[:tfp_sh], exo[:tfp_sh]]       = m[:σ_tfp]^2
    QQ[exo[:gdpdef_sh], exo[:gdpdef_sh]] = m[:σ_gdpdef]^2
    QQ[exo[:pce_sh], exo[:pce_sh]]       = m[:σ_pce]^2

    # These lines set the standard deviations for the anticipated shocks. They
    # are here no longer calibrated to the std dev of contemporaneous shocks,
    # as we had in 904
    if shocks
        for i = 1:anticipated_shocks(m)
            ZZ[obs[symbol("R_n$i")], :] = ZZ[obs[:R_n], :]*(TTT^i)
            DD[obs[symbol("R_n$i")]] = m[:Rstarn]
            QQ[exo[symbol("rm_shl$i")], exo[symbol("rm_shl$i")]] = m[symbol("σ_rm$i")]^2
        end
    end

    return ZZ, DD, QQ, EE, MM
end
