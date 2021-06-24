"""
```
eqcond(m::Model1002)
```

Expresses the equilibrium conditions in canonical form using Γ0, Γ1, C, Ψ, and Π matrices.
Using the mappings of states/equations to integers defined in m1002.jl, coefficients are
specified in their proper positions.

### Outputs

* `Γ0` (`n_states` x `n_states`) holds coefficients of current time states.
* `Γ1` (`n_states` x `n_states`) holds coefficients of lagged states.
* `C`  (`n_states` x `1`) is a vector of constants
* `Ψ`  (`n_states` x `n_shocks_exogenous`) holds coefficients of iid shocks.
* `Π`  (`n_states` x `n_states_expectational`) holds coefficients of expectational states.
"""
function eqcond(m::Model1002)
    return eqcond(m, 1)
end

function eqcond(m::Model1002, reg::Int)
    endo = m.endogenous_states
    exo  = m.exogenous_shocks
    ex   = m.expected_shocks
    eq   = m.equilibrium_conditions

    Γ0 = zeros(n_states(m), n_states(m))
    Γ1 = zeros(n_states(m), n_states(m))
    C  = zeros(n_states(m))
    Ψ  = zeros(n_states(m), n_shocks_exogenous(m))
    Π  = zeros(n_states(m), n_shocks_expectational(m))

    for para in m.parameters
        if !isempty(para.regimes)
            if (haskey(get_settings(m), :model2para_regime) ? haskey(get_setting(m, :model2para_regime), para.key) : false)
                ModelConstructors.toggle_regime!(para, reg, get_setting(m, :model2para_regime)[para.key])
            else
                ModelConstructors.toggle_regime!(para, reg)
            end
        end
    end

    ### ENDOGENOUS STATES ###

    ### 1. Consumption Euler Equation

    # Sticky prices and wages
    Γ0[eq[:eq_euler], endo[:c_t]]  = 1.
    Γ0[eq[:eq_euler], endo[:R_t]]  = (1 - m[:h]*exp(-m[:z_star]))/(m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))
    Γ0[eq[:eq_euler], endo[:b_t]]  = -1.
    Γ0[eq[:eq_euler], endo[:Eπ_t]] = -(1 - m[:h]*exp(-m[:z_star]))/(m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))
    Γ0[eq[:eq_euler], endo[:z_t]]  = (m[:h]*exp(-m[:z_star]))/(1 + m[:h]*exp(-m[:z_star]))
    Γ0[eq[:eq_euler], endo[:Ec_t]] = -1/(1 + m[:h]*exp(-m[:z_star]))
    Γ0[eq[:eq_euler], endo[:Ez_t]] = -1/(1 + m[:h]*exp(-m[:z_star]))
    Γ0[eq[:eq_euler], endo[:L_t]]  = -(m[:σ_c] - 1)*m[:wl_c]/(m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))
    Γ0[eq[:eq_euler], endo[:EL_t]] = (m[:σ_c] - 1)*m[:wl_c]/(m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))
    Γ1[eq[:eq_euler], endo[:c_t]]  = (m[:h]*exp(-m[:z_star]))/(1 + m[:h]*exp(-m[:z_star]))

    if subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63", "ss64", "ss65", "ss66", "ss67", "ss68", "ss69", "ss70", "ss71", "ss72", "ss73", "ss74", "ss75", "ss76", "ss77", "ss78", "ss79", "ss80", "ss81", "ss82", "ss83", "ss84", "ss85", "ss86"]
        Γ0[eq[:eq_euler], endo[:φ_t]]  = -(m[:σ_c] - 1)*m[:wl_c]/(m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))
        Γ0[eq[:eq_euler], endo[:Eφ_t]] = (m[:σ_c] - 1)*m[:wl_c]/(m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))
    end

    # Flexible prices and wages
    Γ0[eq[:eq_euler_f], endo[:c_f_t]]  = 1.
    Γ0[eq[:eq_euler_f], endo[:r_f_t]]  = (1 - m[:h]*exp(-m[:z_star]))/(m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))
    Γ0[eq[:eq_euler_f], endo[:b_t]]    = -1.
    Γ0[eq[:eq_euler_f], endo[:z_t]]    = (m[:h]*exp(-m[:z_star]))/(1 + m[:h]*exp(-m[:z_star]))
    Γ0[eq[:eq_euler_f], endo[:Ec_f_t]] = -1/(1 + m[:h]*exp(-m[:z_star]))
    Γ0[eq[:eq_euler_f], endo[:Ez_t]]   = -1/(1 + m[:h]*exp(-m[:z_star]))
    Γ0[eq[:eq_euler_f], endo[:L_f_t]]  = -(m[:σ_c] - 1)*m[:wl_c]/(m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))
    Γ0[eq[:eq_euler_f], endo[:EL_f_t]] = (m[:σ_c] - 1)*m[:wl_c]/(m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))
    Γ1[eq[:eq_euler_f], endo[:c_f_t]]  = (m[:h]*exp(-m[:z_star]))/(1 + m[:h]*exp(-m[:z_star]))

    if subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63", "ss64", "ss65", "ss66", "ss67", "ss68", "ss69", "ss70", "ss71", "ss72", "ss73", "ss74", "ss75", "ss76", "ss77", "ss78", "ss79", "ss80", "ss81", "ss82", "ss83", "ss84", "ss85", "ss86"]
        Γ0[eq[:eq_euler_f], endo[:φ_t]]  = -(m[:σ_c] - 1)*m[:wl_c]/(m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))
        Γ0[eq[:eq_euler_f], endo[:Eφ_t]] = (m[:σ_c] - 1)*m[:wl_c]/(m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))
    end

    ### 2. Investment Euler Equation

    # Sticky prices and wages
    Γ0[eq[:eq_inv], endo[:qk_t]] = -1/(m[:S′′]*exp(2.0*m[:z_star])*(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star])))
    Γ0[eq[:eq_inv], endo[:i_t]]  = 1.
    Γ0[eq[:eq_inv], endo[:z_t]]  = 1/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
    Γ1[eq[:eq_inv], endo[:i_t]]  = 1/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
    Γ0[eq[:eq_inv], endo[:Ei_t]] = -m[:β]*exp((1 - m[:σ_c])*m[:z_star])/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
    Γ0[eq[:eq_inv], endo[:Ez_t]] = -m[:β]*exp((1 - m[:σ_c])*m[:z_star])/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
    Γ0[eq[:eq_inv], endo[:μ_t]]  = -1.

    # Flexible prices and wages
    Γ0[eq[:eq_inv_f], endo[:qk_f_t]] = -1/(m[:S′′]*exp(2*m[:z_star])*(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star])))
    Γ0[eq[:eq_inv_f], endo[:i_f_t]]  = 1.
    Γ0[eq[:eq_inv_f], endo[:z_t]]    = 1/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
    Γ1[eq[:eq_inv_f], endo[:i_f_t]]  = 1/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
    Γ0[eq[:eq_inv_f], endo[:Ei_f_t]] = -m[:β]*exp((1 - m[:σ_c])*m[:z_star])/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
    Γ0[eq[:eq_inv_f], endo[:Ez_t]]   = -m[:β]*exp((1 - m[:σ_c])*m[:z_star])/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
    Γ0[eq[:eq_inv_f], endo[:μ_t]]    = -1.

    ### 3. Financial Friction Block

    # Return to capital
    # Sticky prices and wages
    Γ0[eq[:eq_capval], endo[:Rktil_t]] = 1.
    Γ0[eq[:eq_capval], endo[:π_t]]      = -1.
    Γ0[eq[:eq_capval], endo[:rk_t]]     = -m[:r_k_star]/(1 + m[:r_k_star] - m[:δ])
    Γ0[eq[:eq_capval], endo[:qk_t]]     = -(1 - m[:δ])/(1 + m[:r_k_star] - m[:δ])
    Γ1[eq[:eq_capval], endo[:qk_t]]     = -1.

    # Spreads
    # Sticky prices and wages
    Γ0[eq[:eq_spread], endo[:ERktil_t]] = 1.
    Γ0[eq[:eq_spread], endo[:R_t]]       = -1.
    Γ0[eq[:eq_spread], endo[:b_t]]       = (m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))/(1 - m[:h]*exp(-m[:z_star]))
    Γ0[eq[:eq_spread], endo[:qk_t]]      = -m[:ζ_spb]
    Γ0[eq[:eq_spread], endo[:kbar_t]]    = -m[:ζ_spb]
    Γ0[eq[:eq_spread], endo[:n_t]]       = m[:ζ_spb]
    Γ0[eq[:eq_spread], endo[:σ_ω_t]]     = -1.
    Γ0[eq[:eq_spread], endo[:μ_e_t]]     = -1.

    # Flexible prices and wages
    Γ0[eq[:eq_spread_f], endo[:ERktil_f_t]] = 1.
    Γ0[eq[:eq_spread_f], endo[:r_f_t]]       = -1.
    Γ0[eq[:eq_spread_f], endo[:b_t]]       = (m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))/(1 - m[:h]*exp(-m[:z_star]))
    Γ0[eq[:eq_spread_f], endo[:qk_f_t]]      = -m[:ζ_spb]
    Γ0[eq[:eq_spread_f], endo[:kbar_f_t]]    = -m[:ζ_spb]
    Γ0[eq[:eq_spread_f], endo[:n_f_t]]       = m[:ζ_spb]
    Γ0[eq[:eq_spread_f], endo[:σ_ω_t]]     = -1.
    Γ0[eq[:eq_spread_f], endo[:μ_e_t]]     = -1.

    # n evol
    # Sticky prices and wages
    Γ0[eq[:eq_nevol], endo[:n_t]]      = 1.
    Γ0[eq[:eq_nevol], endo[:γ_t]]      = -1.
    Γ0[eq[:eq_nevol], endo[:z_t]]      = m[:γ_star]*m[:vstar]/m[:nstar]
    Γ0[eq[:eq_nevol], endo[:Rktil_t]] = -m[:ζ_nRk]
    Γ0[eq[:eq_nevol], endo[:π_t]]      = (m[:ζ_nRk] - m[:ζ_nR])
    Γ1[eq[:eq_nevol], endo[:σ_ω_t]]    = -m[:ζ_nσ_ω]/m[:ζ_spσ_ω]
    Γ1[eq[:eq_nevol], endo[:μ_e_t]]    = -m[:ζ_nμ_e]/m[:ζ_spμ_e]
    Γ1[eq[:eq_nevol], endo[:qk_t]]     = m[:ζ_nqk]
    Γ1[eq[:eq_nevol], endo[:kbar_t]]   = m[:ζ_nqk]
    Γ1[eq[:eq_nevol], endo[:n_t]]      = m[:ζ_nn]
    Γ1[eq[:eq_nevol], endo[:R_t]]      = -m[:ζ_nR]
    Γ1[eq[:eq_nevol], endo[:b_t]]      = m[:ζ_nR]*((m[:σ_c]*(1.0+m[:h]*exp(-m[:z_star])))/(1.0-m[:h]*exp(-m[:z_star])))

    # Flexible prices and wages
    Γ0[eq[:eq_nevol_f], endo[:n_f_t]]      = 1.
    Γ0[eq[:eq_nevol_f], endo[:z_t]]      = m[:γ_star]*m[:vstar]/m[:nstar]
    Γ0[eq[:eq_nevol_f], endo[:Rktil_f_t]] = -m[:ζ_nRk]
    Γ1[eq[:eq_nevol_f], endo[:σ_ω_t]]    = -m[:ζ_nσ_ω]/m[:ζ_spσ_ω]
    Γ1[eq[:eq_nevol_f], endo[:μ_e_t]]    = -m[:ζ_nμ_e]/m[:ζ_spμ_e]
    Γ1[eq[:eq_nevol_f], endo[:qk_f_t]]     = m[:ζ_nqk]
    Γ1[eq[:eq_nevol_f], endo[:kbar_f_t]]   = m[:ζ_nqk]
    Γ1[eq[:eq_nevol_f], endo[:n_f_t]]      = m[:ζ_nn]
    Γ1[eq[:eq_nevol_f], endo[:r_f_t]]      = -m[:ζ_nR]
    Γ1[eq[:eq_nevol_f], endo[:b_t]]      = m[:ζ_nR]*((m[:σ_c]*(1.0+m[:h]*exp(-m[:z_star])))/(1.0-m[:h]*exp(-m[:z_star])))

    # Flexible prices and wages - ASSUME NO FINANCIAL FRICTIONS
    Γ0[eq[:eq_capval_f], endo[:Rktil_f_t]] = 1.
    Γ0[eq[:eq_capval_f], endo[:rk_f_t]]     = -m[:r_k_star]/(m[:r_k_star]+1-m[:δ])
    Γ0[eq[:eq_capval_f], endo[:qk_f_t]]     = -(1-m[:δ])/(m[:r_k_star]+1-m[:δ])
    Γ1[eq[:eq_capval_f], endo[:qk_f_t]]     = -1.

    ### 4. Aggregate Production Function

    # Sticky prices and wages
    Γ0[eq[:eq_output], endo[:y_t]] =  1.
    Γ0[eq[:eq_output], endo[:k_t]] = -m[:Φ]*m[:α]
    Γ0[eq[:eq_output], endo[:L_t]] = -m[:Φ]*(1 - m[:α])

    # Flexible prices and wages
    Γ0[eq[:eq_output_f], endo[:y_f_t]] =  1.
    Γ0[eq[:eq_output_f], endo[:k_f_t]] = -m[:Φ]*m[:α]
    Γ0[eq[:eq_output_f], endo[:L_f_t]] = -m[:Φ]*(1 - m[:α])

    ### 5. Capital Utilization

    # Sticky prices and wages
    Γ0[eq[:eq_caputl], endo[:k_t]]    =  1.
    Γ1[eq[:eq_caputl], endo[:kbar_t]] =  1.
    Γ0[eq[:eq_caputl], endo[:z_t]]    = 1.
    Γ0[eq[:eq_caputl], endo[:u_t]]    = -1.

    # Flexible prices and wages
    Γ0[eq[:eq_caputl_f], endo[:k_f_t]]    =  1.
    Γ1[eq[:eq_caputl_f], endo[:kbar_f_t]] =  1.
    Γ0[eq[:eq_caputl_f], endo[:z_t]]      = 1.
    Γ0[eq[:eq_caputl_f], endo[:u_f_t]]    = -1.

    ### 6. Rental Rate of Capital

    # Sticky prices and wages
    Γ0[eq[:eq_capsrv], endo[:u_t]]  = 1.
    Γ0[eq[:eq_capsrv], endo[:rk_t]] = -(1 - m[:ppsi])/m[:ppsi]

    # Flexible prices and wages
    Γ0[eq[:eq_capsrv_f], endo[:u_f_t]]  = 1.
    Γ0[eq[:eq_capsrv_f], endo[:rk_f_t]] = -(1 - m[:ppsi])/m[:ppsi]

    ### 7. Evolution of Capital

    # Sticky prices and wages
    Γ0[eq[:eq_capev], endo[:kbar_t]] = 1.
    Γ1[eq[:eq_capev], endo[:kbar_t]] = 1 - m[:istar]/m[:kbarstar]
    Γ0[eq[:eq_capev], endo[:z_t]]    = 1 - m[:istar]/m[:kbarstar]
    Γ0[eq[:eq_capev], endo[:i_t]]    = -m[:istar]/m[:kbarstar]
    Γ0[eq[:eq_capev], endo[:μ_t]]    = -m[:istar]*m[:S′′]*exp(2*m[:z_star])*(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))/m[:kbarstar]

    # Flexible prices and wages
    Γ0[eq[:eq_capev_f], endo[:kbar_f_t]] = 1.
    Γ1[eq[:eq_capev_f], endo[:kbar_f_t]] = 1 - m[:istar]/m[:kbarstar]
    Γ0[eq[:eq_capev_f], endo[:z_t]]      = 1 - m[:istar]/m[:kbarstar]
    Γ0[eq[:eq_capev_f], endo[:i_f_t]]    = -m[:istar]/m[:kbarstar]
    Γ0[eq[:eq_capev_f], endo[:μ_t]]      = -m[:istar]*m[:S′′]*exp(2*m[:z_star])*(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))/m[:kbarstar]

    ### 8. Price Markup

    # Sticky prices and wages
    Γ0[eq[:eq_mkupp], endo[:mc_t]] =  1.
    Γ0[eq[:eq_mkupp], endo[:w_t]]  = -1.
    Γ0[eq[:eq_mkupp], endo[:L_t]]  = -m[:α]
    Γ0[eq[:eq_mkupp], endo[:k_t]]  =  m[:α]

    # Flexible prices and wages
    Γ0[eq[:eq_mkupp_f], endo[:w_f_t]] = 1.
    Γ0[eq[:eq_mkupp_f], endo[:L_f_t]] =  m[:α]
    Γ0[eq[:eq_mkupp_f], endo[:k_f_t]] =  -m[:α]

    ### 9. Phillips Curve

    # Sticky prices and wages
    Γ0[eq[:eq_phlps], endo[:π_t]]  = 1.
    Γ0[eq[:eq_phlps], endo[:mc_t]] =  -((1 - m[:ζ_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))*
        (1 - m[:ζ_p]))/(m[:ζ_p]*((m[:Φ]- 1)*m[:ϵ_p] + 1))/(1 + m[:ι_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
    Γ1[eq[:eq_phlps], endo[:π_t]]  = m[:ι_p]/(1 + m[:ι_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
    Γ0[eq[:eq_phlps], endo[:Eπ_t]] = -m[:β]*exp((1 - m[:σ_c])*m[:z_star])/(1 + m[:ι_p]*m[:β]*
        exp((1 - m[:σ_c])*m[:z_star]))

    # Comment out for counterfactual with no price mark up shock
    Γ0[eq[:eq_phlps], endo[:λ_f_t]] = -1.

    # Flexible prices and wages not necessary

    ### 10. Rental Rate of Capital

    # Sticky prices and wages
    Γ0[eq[:eq_caprnt], endo[:rk_t]] = 1.
    Γ0[eq[:eq_caprnt], endo[:k_t]]  = 1.
    Γ0[eq[:eq_caprnt], endo[:L_t]]  = -1.
    Γ0[eq[:eq_caprnt], endo[:w_t]]  = -1.

    # Flexible prices and wages
    Γ0[eq[:eq_caprnt_f], endo[:rk_f_t]] = 1.
    Γ0[eq[:eq_caprnt_f], endo[:k_f_t]]  = 1.
    Γ0[eq[:eq_caprnt_f], endo[:L_f_t]]  = -1.
    Γ0[eq[:eq_caprnt_f], endo[:w_f_t]]  = -1.

    ### 11. Marginal Substitution

    # Sticky prices and wages
    Γ0[eq[:eq_msub], endo[:μ_ω_t]] = 1.
    Γ0[eq[:eq_msub], endo[:L_t]]   = m[:ν_l]
    Γ0[eq[:eq_msub], endo[:c_t]]   = 1/(1 - m[:h]*exp(-m[:z_star]))
    Γ1[eq[:eq_msub], endo[:c_t]]   = m[:h]*exp(-m[:z_star])/(1 - m[:h]*exp(-m[:z_star]))
    Γ0[eq[:eq_msub], endo[:z_t]]   = m[:h]*exp(-m[:z_star]) /(1 - m[:h]*exp(-m[:z_star]))
    Γ0[eq[:eq_msub], endo[:w_t]]   = -1.

    if subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63", "ss64", "ss65", "ss66", "ss67", "ss68", "ss69", "ss70", "ss71", "ss72", "ss73", "ss74", "ss75", "ss76", "ss77", "ss78", "ss79", "ss80", "ss81", "ss82", "ss83", "ss84", "ss85", "ss86"]
        Γ0[eq[:eq_msub], endo[:φ_t]] = m[:ν_l]
    end

    # Flexible prices and wages
    Γ0[eq[:eq_msub_f], endo[:w_f_t]] = -1.
    Γ0[eq[:eq_msub_f], endo[:L_f_t]] = m[:ν_l]
    Γ0[eq[:eq_msub_f], endo[:c_f_t]] = 1/(1 - m[:h]*exp(-m[:z_star]))
    Γ1[eq[:eq_msub_f], endo[:c_f_t]] = m[:h]*exp(-m[:z_star])/(1 - m[:h]*exp(-m[:z_star]))
    Γ0[eq[:eq_msub_f], endo[:z_t]]   = m[:h]*exp(-m[:z_star])/(1 - m[:h]*exp(-m[:z_star]))

    if subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63", "ss64", "ss65", "ss66", "ss67", "ss68", "ss69", "ss70", "ss71", "ss72", "ss73", "ss74", "ss75", "ss76", "ss77", "ss78", "ss79", "ss80", "ss81", "ss82", "ss83", "ss84", "ss85", "ss86"]
        Γ0[eq[:eq_msub_f], endo[:φ_t]] = m[:ν_l]
    end

    ### 12. Evolution of Wages

    # Sticky prices and wages
    Γ0[eq[:eq_wage], endo[:w_t]]   = 1
    Γ0[eq[:eq_wage], endo[:μ_ω_t]] = (1 - m[:ζ_w]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))*
        (1 - m[:ζ_w])/(m[:ζ_w]*((m[:λ_w] - 1)*m[:ϵ_w] + 1))/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
    Γ0[eq[:eq_wage], endo[:π_t]]   = (1 + m[:ι_w]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
    Γ1[eq[:eq_wage], endo[:w_t]]   = 1/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
    Γ0[eq[:eq_wage], endo[:z_t]]   = 1/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
    Γ1[eq[:eq_wage], endo[:π_t]]   = m[:ι_w]/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
    Γ0[eq[:eq_wage], endo[:Ew_t]]  = -m[:β]*exp((1 - m[:σ_c])*m[:z_star])/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
    Γ0[eq[:eq_wage], endo[:Ez_t]]  = -m[:β]*exp((1 - m[:σ_c])*m[:z_star])/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
    Γ0[eq[:eq_wage], endo[:Eπ_t]]  = -m[:β]*exp((1 - m[:σ_c])*m[:z_star])/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
    Γ0[eq[:eq_wage], endo[:λ_w_t]] = -1.

    # Flexible prices and wages not necessary

    ### 13. Monetary Policy Rule
    # Sticky prices and wages
    Γ0[eq[:eq_mp], endo[:R_t]]      = 1.
    Γ1[eq[:eq_mp], endo[:R_t]]      = m[:ρ]
    Γ0[eq[:eq_mp], endo[:π_t]]      = -(1 - m[:ρ])*m[:ψ1]
    Γ0[eq[:eq_mp], endo[:π_star_t]] = (1 - m[:ρ])*m[:ψ1]
    Γ0[eq[:eq_mp], endo[:y_t]]      = -(1 - m[:ρ])*m[:ψ2] - m[:ψ3]
    Γ0[eq[:eq_mp], endo[:y_f_t]]    = (1 - m[:ρ])*m[:ψ2] + m[:ψ3]
    Γ1[eq[:eq_mp], endo[:y_t]]      = -m[:ψ3]
    Γ1[eq[:eq_mp], endo[:y_f_t]]    = m[:ψ3]
    Γ0[eq[:eq_mp], endo[:rm_t]]     = -1.

    # Flexible prices and wages not necessary

    ### 14. Resource Constraint

    # Sticky prices and wages
    Γ0[eq[:eq_res], endo[:y_t]] = 1.
    Γ0[eq[:eq_res], endo[:g_t]] = -m[:g_star]
    Γ0[eq[:eq_res], endo[:c_t]] = -m[:cstar]/m[:ystar]
    Γ0[eq[:eq_res], endo[:i_t]] = -m[:istar]/m[:ystar]
    Γ0[eq[:eq_res], endo[:u_t]] = -m[:r_k_star]*m[:kstar]/m[:ystar]

    # Flexible prices and wages
    Γ0[eq[:eq_res_f], endo[:y_f_t]] = 1.
    Γ0[eq[:eq_res_f], endo[:g_t]]   = -m[:g_star]
    Γ0[eq[:eq_res_f], endo[:c_f_t]] = -m[:cstar]/m[:ystar]
    Γ0[eq[:eq_res_f], endo[:i_f_t]] = -m[:istar]/m[:ystar]
    Γ0[eq[:eq_res_f], endo[:u_f_t]] = -m[:r_k_star]*m[:kstar]/m[:ystar]

    ### 15. Extra States
    # These aren't strictly necessary, but they track lags or simplify the equations

    # π_t1
    Γ0[eq[:eq_π1], endo[:π_t1]] = 1.
    Γ1[eq[:eq_π1], endo[:π_t]]  = 1.

    # π_t2
    Γ0[eq[:eq_π2], endo[:π_t2]] = 1.
    Γ1[eq[:eq_π2], endo[:π_t1]] = 1.

    # π_a
    Γ0[eq[:eq_π_a], endo[:π_a_t]] = 1.
    Γ0[eq[:eq_π_a], endo[:π_t]]   = -1.
    Γ0[eq[:eq_π_a], endo[:π_t1]]  = -1.
    Γ0[eq[:eq_π_a], endo[:π_t2]]  = -1.
    Γ1[eq[:eq_π_a], endo[:π_t2]]  = 1.

    # Rt1
    Γ0[eq[:eq_Rt1], endo[:R_t1]] = 1.
    Γ1[eq[:eq_Rt1], endo[:R_t]]  = 1.

    # Ez_t
    Γ0[eq[:eq_Ez], endo[:Ez_t]]   = 1.
    Γ0[eq[:eq_Ez], endo[:ztil_t]] = -(m[:ρ_ztil]-1)/(1-m[:α])
    Γ0[eq[:eq_Ez], endo[:zp_t]]   = -m[:ρ_z_p]

    # Eφ_t
    if subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63", "ss64", "ss65", "ss66", "ss67", "ss68", "ss69", "ss70", "ss71", "ss72", "ss73", "ss74", "ss75", "ss76", "ss77", "ss78", "ss79", "ss80", "ss81", "ss82", "ss83", "ss84", "ss85", "ss86"]
        Γ0[eq[:eq_Eφ], endo[:Eφ_t]] = 1.
        Γ0[eq[:eq_Eφ], endo[:φ_t]]  = -m[:ρ_φ]
    end

    ### EXOGENOUS SHOCKS ###

    # Neutral technology
    Γ0[eq[:eq_z], endo[:z_t]]    = 1.
    Γ1[eq[:eq_z], endo[:ztil_t]] = (m[:ρ_ztil] - 1)/(1 - m[:α])
    Γ0[eq[:eq_z], endo[:zp_t]]   = -1.
    Ψ[eq[:eq_z], exo[:ztil_sh]]     = 1/(1 - m[:α])

    Γ0[eq[:eq_ztil], endo[:ztil_t]] = 1.
    Γ1[eq[:eq_ztil], endo[:ztil_t]] = m[:ρ_ztil]
    Ψ[eq[:eq_ztil], exo[:ztil_sh]]     = 1.

    if subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63", "ss64", "ss65", "ss66", "ss67", "ss68", "ss69", "ss70", "ss71", "ss72", "ss73", "ss74", "ss75", "ss76", "ss77", "ss78", "ss79", "ss80", "ss81", "ss82", "ss83", "ss84", "ss85", "ss86"]
        # Ez_t
        Γ0[eq[:eq_Ez], endo[:ziid_t]] = -(m[:ρ_ziid]-1)/(1-m[:α])

        # Neutral technology
        Γ0[eq[:eq_z], endo[:ziid_t]]  = -1 / (1 - m[:α])
        Γ1[eq[:eq_z], endo[:ziid_t]]  = -1 / (1 - m[:α])

        # AR(1) for ziid
        Γ0[eq[:eq_ziid], endo[:ziid_t]] = 1.
        Γ1[eq[:eq_ziid], endo[:ziid_t]] = m[:ρ_ziid]
        Ψ[eq[:eq_ziid], exo[:ziid_sh]]  = 1.
    end

    # Long-run changes to productivity
    Γ0[eq[:eq_zp], endo[:zp_t]] = 1.
    Γ1[eq[:eq_zp], endo[:zp_t]] = m[:ρ_z_p]
    Ψ[eq[:eq_zp], exo[:zp_sh]]  = 1.

    # Government spending
    Γ0[eq[:eq_g], endo[:g_t]] = 1.
    Γ1[eq[:eq_g], endo[:g_t]] = m[:ρ_g]
    Ψ[eq[:eq_g], exo[:g_sh]]  = 1.
    Ψ[eq[:eq_g], exo[:ztil_sh]]  = m[:η_gz]

    # Asset shock
    Γ0[eq[:eq_b], endo[:b_t]] = 1.
    Γ1[eq[:eq_b], endo[:b_t]] = m[:ρ_b]
    Ψ[eq[:eq_b], exo[:b_sh]]  = 1.

    if subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63", "ss64", "ss65", "ss66", "ss67", "ss68", "ss69", "ss70", "ss71", "ss72", "ss73", "ss74", "ss75", "ss76", "ss77", "ss78", "ss79", "ss80", "ss81", "ss82", "ss83", "ss84", "ss85", "ss86"]
        # iid shock to Euler equation
        Γ0[eq[:eq_biidc], endo[:biidc_t]] = 1.
        Γ1[eq[:eq_biidc], endo[:biidc_t]] = m[:ρ_biidc] # c b/c will only affect consumption
        Ψ[eq[:eq_biidc], exo[:biidc_sh]]  = 1.

        Γ0[eq[:eq_euler], endo[:biidc_t]]   = -1.
        Γ0[eq[:eq_euler_f], endo[:biidc_t]] = -1.
    end

    # Investment-specific technology
    Γ0[eq[:eq_μ], endo[:μ_t]] = 1.
    Γ1[eq[:eq_μ], endo[:μ_t]] = m[:ρ_μ]
    Ψ[eq[:eq_μ], exo[:μ_sh]]  = 1.

    # Price mark-up shock
    Γ0[eq[:eq_λ_f], endo[:λ_f_t]]  = 1.
    Γ1[eq[:eq_λ_f], endo[:λ_f_t]]  = m[:ρ_λ_f]
    Γ1[eq[:eq_λ_f], endo[:λ_f_t1]] = -m[:η_λ_f]
    Ψ[eq[:eq_λ_f], exo[:λ_f_sh]]   = 1.

    Γ0[eq[:eq_λ_f1], endo[:λ_f_t1]] = 1.
    Ψ[eq[:eq_λ_f1], exo[:λ_f_sh]]   = 1.

    if subspec(m) in ["ss86"]
        Ψ[eq[:eq_λ_f], exo[:λ_f_iid_sh]] = 1.0
    end

    # Wage mark-up shock
    Γ0[eq[:eq_λ_w], endo[:λ_w_t]]  = 1.
    Γ1[eq[:eq_λ_w], endo[:λ_w_t]]  = m[:ρ_λ_w]
    Γ1[eq[:eq_λ_w], endo[:λ_w_t1]] = -m[:η_λ_w]
    Ψ[eq[:eq_λ_w], exo[:λ_w_sh]]   = 1.

    Γ0[eq[:eq_λ_w1], endo[:λ_w_t1]] = 1.
    Ψ[eq[:eq_λ_w1], exo[:λ_w_sh]]   = 1.

    # Monetary policy shock
    noant = haskey(m.settings, :remove_rm_t_shocks) &&
            reg >= get_setting(m, :remove_rm_t_shocks) ? 0.0 : 1.0

    Γ0[eq[:eq_rm], endo[:rm_t]] = 1.
    Γ1[eq[:eq_rm], endo[:rm_t]] = noant * m[:ρ_rm]
    Ψ[eq[:eq_rm], exo[:rm_sh]]  = noant

    # Labor preference shock
    if subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63", "ss64", "ss65", "ss66", "ss67", "ss68", "ss69", "ss70", "ss71", "ss72", "ss73", "ss74", "ss75", "ss76", "ss77", "ss78", "ss79", "ss80", "ss81", "ss82", "ss83", "ss84", "ss85", "ss86"]
        # Eφ_t
        Γ0[eq[:eq_Eφ], endo[:φ_t]] = -m[:ρ_φ]

        # AR(1) process
        Γ0[eq[:eq_φ], endo[:φ_t]] = 1.
        Γ1[eq[:eq_φ], endo[:φ_t]] = m[:ρ_φ]
        Ψ[eq[:eq_φ], exo[:φ_sh]]  = 1.
    end

    # COVID counterparts for standard business cycle shocks
    if subspec(m) in ["ss67", "ss68", "ss69", "ss70", "ss71", "ss72", "ss73", "ss74", "ss75", "ss76", "ss77", "ss78", "ss80", "ss82", "ss83"]
        # TODO: check if there are any expectational terms to which
        # we need to account for (see Eφ_t)
        Γ0[eq[:eq_g_covid], endo[:g_covid_t]] = 1.
        Γ1[eq[:eq_g_covid], endo[:g_covid_t]] = m[:ρ_g_covid]
        Ψ[eq[:eq_g_covid], exo[:g_covid_sh]]  = 1.
        Γ0[eq[:eq_g], endo[:g_covid_t]]       = -1.

        Γ0[eq[:eq_μ_covid], endo[:μ_covid_t]] = 1.
        Γ1[eq[:eq_μ_covid], endo[:μ_covid_t]] = m[:ρ_μ_covid]
        Ψ[eq[:eq_μ_covid],  exo[:μ_covid_sh]] = 1.
        Γ0[eq[:eq_μ], endo[:μ_covid_t]]       = -1.

        Γ0[eq[:eq_λ_f_covid], endo[:λ_f_covid_t]] = 1.
        Γ1[eq[:eq_λ_f_covid], endo[:λ_f_covid_t]] = m[:ρ_λ_f_covid]
        Ψ[eq[:eq_λ_f_covid], exo[:λ_f_covid_sh]]  = 1.
        Γ0[eq[:eq_λ_f], endo[:λ_f_covid_t]]       = -1.

        Γ0[eq[:eq_σ_ω_covid], endo[:σ_ω_covid_t]] = 1.
        Γ1[eq[:eq_σ_ω_covid], endo[:σ_ω_covid_t]] = m[:ρ_σ_w_covid]
        Ψ[eq[:eq_σ_ω_covid], exo[:σ_ω_covid_sh]]  = 1.
        Γ0[eq[:eq_σ_ω], endo[:σ_ω_covid_t]]       = -1.
    end

    if subspec(m) in ["ss69", "ss70", "ss73", "ss74", "ss77", "ss78"]
        Γ0[eq[:eq_zp_covid], endo[:zp_covid_t]] = 1.
        Γ1[eq[:eq_zp_covid], endo[:zp_covid_t]] = m[:ρ_z_p_covid]
        Ψ[eq[:eq_zp_covid], exo[:zp_covid_sh]]  = 1.
        Γ0[eq[:eq_zp], endo[:zp_covid_t]]       = -1.
    end

    ### Financial frictions

    # Standard deviation of capital shock to entrepreneurs
    Γ0[eq[:eq_σ_ω], endo[:σ_ω_t]] = 1.
    Γ1[eq[:eq_σ_ω], endo[:σ_ω_t]] = m[:ρ_σ_w]
    Ψ[eq[:eq_σ_ω], exo[:σ_ω_sh]]  = 1.

    # Exogenous bankruptcy costs
    Γ0[eq[:eq_μ_e], endo[:μ_e_t]] = 1.
    Γ1[eq[:eq_μ_e], endo[:μ_e_t]] = m[:ρ_μ_e]
    Ψ[eq[:eq_μ_e], exo[:μ_e_sh]]  = 1.

    # Fraction of entrepreneurs surviving period t
    Γ0[eq[:eq_γ], endo[:γ_t]] = 1.
    Γ1[eq[:eq_γ], endo[:γ_t]] = m[:ρ_γ]
    Ψ[eq[:eq_γ], exo[:γ_sh]]  = 1.

    # Long-term inflation expectations
    Γ0[eq[:eq_π_star], endo[:π_star_t]] = 1.
    Γ1[eq[:eq_π_star], endo[:π_star_t]] = m[:ρ_π_star]
    Ψ[eq[:eq_π_star], exo[:π_star_sh]]  = 1.

    # Anticipated policy shocks
    if n_mon_anticipated_shocks(m) > 0

        noant = haskey(m.settings, :remove_rm_shocks) &&
            reg >= get_setting(m, :remove_rm_shocks) ? 0.0 : 1.0

        # This section adds the anticipated shocks. There is one state for all the
        # anticipated shocks that will hit in a given period (i.e. rm_tl2 holds those that
        # will hit in two periods), and the equations are set up so that rm_tl2 last period
        # will feed into rm_tl1 this period (and so on for other numbers), and last period's
        # rm_tl1 will feed into the rm_t process (and affect the Taylor Rule this period).
        Γ1[eq[:eq_rm], endo[:rm_tl1]]   = noant
        Γ0[eq[:eq_rml1], endo[:rm_tl1]] = 1.
        Ψ[eq[:eq_rml1], exo[:rm_shl1]]  = noant

        if n_mon_anticipated_shocks(m) > 1
            for i = 2:n_mon_anticipated_shocks(m)
                Γ1[eq[Symbol("eq_rml$(i-1)")], endo[Symbol("rm_tl$i")]] = noant
                Γ0[eq[Symbol("eq_rml$i")], endo[Symbol("rm_tl$i")]]     = 1.
                Ψ[eq[Symbol("eq_rml$i")], exo[Symbol("rm_shl$i")]]      = noant
            end

            #=if (haskey(m.settings, :flexible_ait_policy_change) ? get_setting(m, :flexible_ait_policy_change) : false)
                if get_setting(m, :regime_dates)[reg] >= get_setting(m, :flexible_ait_policy_change_date)
                    Γ1[eq[:eq_rml1], endo[:rm_tl2]] = 0.
                    Γ0[eq[:eq_rml2], endo[:rm_tl2]] = 1.
                    Ψ[eq[:eq_rml2],  exo[:rm_shl2]] = 1.
                end
            end=#
        end
    end

    for (key, val) in get_setting(m, :antshocks)
        ant_eq_mapping   = get_setting(m, :ant_eq_mapping)         # maps antshock key to state variable name
        ant_eq_E_mapping = get_setting(m, :ant_eq_E_mapping)       # maps antshock key to state variable name
        ant_proportion   = get_setting(m, :proportional_antshocks) # is the antshock proportional to a contemporaneous shock?
        ant_contemp_prop = get_setting(m, :contemporaneous_and_proportional_antshocks) # Both proportional and contemporaneous shocks?
        if val > 0
            # This section adds the anticipated shocks. There is one state for all the
            # anticipated shocks that will hit in a given period (i.e. rm_tl2 holds those that
            # will hit in two periods), and the equations are set up so that rm_tl2 last period
            # will feed into rm_tl1 this period (and so on for other numbers), and last period's
            # rm_tl1 will feed into the rm_t process (and affect the Taylor Rule this period).

            Γ1[eq[Symbol("eq_", ant_eq_mapping[key])], endo[Symbol(key, "_tl1")]]   = 1.
            Γ0[eq[Symbol("eq_", key, "l1")], endo[Symbol(key, "_tl1")]] = 1.
            Ψ[eq[Symbol("eq_", key, "l1")], exo[Symbol(key, "_shl1")]]  = 1.

            if key in ant_contemp_prop
                Ψ[eq[Symbol("eq_", key, "l1")], exo[Symbol(key, "_sh")]]   = m[Symbol(:σ_, key, :_prop)]
            elseif key in ant_proportion
                Ψ[eq[Symbol("eq_", key, "l1")], exo[Symbol(key, "_shl1")]] = 0.
                Ψ[eq[Symbol("eq_", key, "l1")], exo[Symbol(key, "_sh")]]   = m[Symbol(:σ_, key, :_prop)]
            end

            if val > 1
                for i = 2:val
                    Γ1[eq[Symbol("eq_", key, "l$(i-1)")], endo[Symbol(key, "_tl$i")]] = 1.
                    Γ0[eq[Symbol("eq_", key, "l$i")], endo[Symbol(key, "_tl$i")]]     = 1.
                    Ψ[eq[Symbol("eq_", key, "l$i")], exo[Symbol(key, "_shl$i")]]      = 1.

                    if key in ant_contemp_prop
                        Ψ[eq[Symbol("eq_", key, "l$i")], exo[Symbol(key, "_sh")]]    = m[Symbol(:σ_, key, :_prop, i)]
                    elseif key in ant_proportion
                        Ψ[eq[Symbol("eq_", key, "l$i")], exo[Symbol(key, "_shl$i")]] = 0.
                        Ψ[eq[Symbol("eq_", key, "l$i")], exo[Symbol(key, "_sh")]]    = m[Symbol(:σ_, key, :_prop, i)]
                    end
                end
            end

            if key in [:z, :ziid] # Handle separately b/c coefficient is not -1.0
                # Ez_t
                Γ0[eq[:eq_Ez], endo[Symbol(key, "_tl1")]] = -1 / (1 - m[:α]) # note z_tl1 = sum of all shocks that will hit next period.
                # so this is the only required line
            elseif haskey(ant_eq_E_mapping, key) # Account for other expected anticipated shocks
                Γ0[eq[Symbol("eq_", ant_eq_E_mapping[key])], endo[Symbol(key, "_tl1")]] = -1.
            end
        end
    end

    ### EXPECTATION ERRORS ###

    ### E(c)

    # Sticky prices and wages
    Γ0[eq[:eq_Ec], endo[:c_t]]  = 1.
    Γ1[eq[:eq_Ec], endo[:Ec_t]] = 1.
    Π[eq[:eq_Ec], ex[:Ec_sh]]   = 1.

    # Flexible prices and wages
    Γ0[eq[:eq_Ec_f], endo[:c_f_t]]  = 1.
    Γ1[eq[:eq_Ec_f], endo[:Ec_f_t]] = 1.
    Π[eq[:eq_Ec_f], ex[:Ec_f_sh]]   = 1.

    ### E(q)

    # Sticky prices and wages
    Γ0[eq[:eq_Eqk], endo[:qk_t]]  = 1.
    Γ1[eq[:eq_Eqk], endo[:Eqk_t]] = 1.
    Π[eq[:eq_Eqk], ex[:Eqk_sh]]   = 1.

    # Flexible prices and wages
    Γ0[eq[:eq_Eqk_f], endo[:qk_f_t]]  = 1.
    Γ1[eq[:eq_Eqk_f], endo[:Eqk_f_t]] = 1.
    Π[eq[:eq_Eqk_f], ex[:Eqk_f_sh]]   = 1.

    ### E(i)

    # Sticky prices and wages
    Γ0[eq[:eq_Ei], endo[:i_t]]  = 1.
    Γ1[eq[:eq_Ei], endo[:Ei_t]] = 1.
    Π[eq[:eq_Ei], ex[:Ei_sh]]   = 1.

    # Flexible prices and wages
    Γ0[eq[:eq_Ei_f], endo[:i_f_t]]  = 1.
    Γ1[eq[:eq_Ei_f], endo[:Ei_f_t]] = 1.
    Π[eq[:eq_Ei_f], ex[:Ei_f_sh]]   = 1.

    ### E(π)

    # Sticky prices and wages
    Γ0[eq[:eq_Eπ], endo[:π_t]]  = 1.
    Γ1[eq[:eq_Eπ], endo[:Eπ_t]] = 1.
    Π[eq[:eq_Eπ], ex[:Eπ_sh]]   = 1.

    ### E(l)

    # Sticky prices and wages
    Γ0[eq[:eq_EL], endo[:L_t]]  = 1.
    Γ1[eq[:eq_EL], endo[:EL_t]] = 1.
    Π[eq[:eq_EL], ex[:EL_sh]]   = 1.

    # Flexible prices and wages
    Γ0[eq[:eq_EL_f], endo[:L_f_t]]  = 1.
    Γ1[eq[:eq_EL_f], endo[:EL_f_t]] = 1.
    Π[eq[:eq_EL_f], ex[:EL_f_sh]]   = 1.

    ### E(rk)

    # Sticky prices and wages
    Γ0[eq[:eq_Erk], endo[:rk_t]]  = 1.
    Γ1[eq[:eq_Erk], endo[:Erk_t]] = 1.
    Π[eq[:eq_Erk], ex[:Erk_sh]]   = 1.

    # Flexible prices and wages
    Γ0[eq[:eq_ERktil_f], endo[:Rktil_f_t]]  = 1.
    Γ1[eq[:eq_ERktil_f], endo[:ERktil_f_t]] = 1.
    Π[eq[:eq_ERktil_f], ex[:ERktil_f_sh]]   = 1.

    ### E(w)

    # Sticky prices and wages
    Γ0[eq[:eq_Ew], endo[:w_t]]  = 1.
    Γ1[eq[:eq_Ew], endo[:Ew_t]] = 1.
    Π[eq[:eq_Ew], ex[:Ew_sh]]   = 1.

    ### E(Rktil)

    # Sticky prices and wages
    Γ0[eq[:eq_ERktil], endo[:Rktil_t]]  = 1.
    Γ1[eq[:eq_ERktil], endo[:ERktil_t]] = 1.
    Π[eq[:eq_ERktil], ex[:ERktil_sh]]    = 1.

   # We additionally need to directly add the equation(s) for pgap, ygap, etc. here rather than just
   # in the altpolicy files b/c
   # (1) Regime-switching won't work otherwise
   # (2) we may want to define the their values at the beginning of the forecast period
   #     rather than just when we start using the alt rule.
   if haskey(m.settings, :add_altpolicy_pgap) ? get_setting(m, :add_altpolicy_pgap) : false
       Γ0[eq[:eq_pgap], endo[:pgap_t]]  =  1.
       if haskey(m.settings, :regime_eqcond_info)
           if reg >= minimum(keys(get_setting(m, :regime_eqcond_info))) &&
               reg <= maximum(keys(get_setting(m, :regime_eqcond_info))) &&
               haskey(m.settings, :pgap_type)
               if get_setting(m, :pgap_type) == :ngdp
                   Γ0[eq[:eq_pgap], endo[:pgap_t]] =  1.
                   Γ0[eq[:eq_pgap], endo[:π_t]]    = -1.
                   Γ1[eq[:eq_pgap], endo[:pgap_t]] =  1.

                   Γ0[eq[:eq_pgap], endo[:y_t]]    = -1.
                   Γ0[eq[:eq_pgap], endo[:z_t]]    = -1.
                   Γ1[eq[:eq_pgap], endo[:y_t]]    = -1.
               elseif get_setting(m, :pgap_type) == :ait
                   Thalf = haskey(get_settings(m), :ait_Thalf) ? get_setting(m, :ait_Thalf) : 10.
                   ρ_ait = exp(log(0.5) / Thalf)
                   Γ0[eq[:eq_pgap], endo[:pgap_t]] = 1.
                   Γ0[eq[:eq_pgap], endo[:π_t]]    = -1.
                   Γ1[eq[:eq_pgap], endo[:pgap_t]] = ρ_ait
               elseif get_setting(m, :pgap_type) == :smooth_ait
                   Thalf = haskey(get_settings(m), :ait_Thalf) ? get_setting(m, :ait_Thalf) : 10.
                   ρ_pgap = exp(log(0.5) / Thalf)
                   Γ0[eq[:eq_pgap], endo[:pgap_t]] = 1.
                   Γ0[eq[:eq_pgap], endo[:π_t]]    = -1.
                   Γ1[eq[:eq_pgap], endo[:pgap_t]] = ρ_pgap
               elseif get_setting(m, :pgap_type) in [:smooth_ait_gdp, :smooth_ait_gdp_alt, :flexible_ait, :rw]
                   Thalf = haskey(get_settings(m), :ait_Thalf) ? get_setting(m, :ait_Thalf) : 10.
                   ρ_pgap = exp(log(0.5) / Thalf)
                   Γ0[eq[:eq_pgap], endo[:pgap_t]] = 1.
                   Γ0[eq[:eq_pgap], endo[:π_t]]    = -1.
                   Γ1[eq[:eq_pgap], endo[:pgap_t]] = ρ_pgap
               end
               if (haskey(get_settings(m), :add_initialize_pgap_ygap_pseudoobs) ?
                   get_setting(m, :add_initialize_pgap_ygap_pseudoobs) : false)
                   Ψ[eq[:eq_pgap], exo[:pgap_sh]] = 1.
               end
           end
       end
       if (haskey(get_settings(m), :add_initialize_pgap_ygap_pseudoobs) ?
                   get_setting(m, :add_initialize_pgap_ygap_pseudoobs) : false)
                   Ψ[eq[:eq_pgap], exo[:pgap_sh]] = 1.
               end
   end

   if haskey(m.settings, :add_altpolicy_ygap) ? get_setting(m, :add_altpolicy_ygap) : false
       Γ0[eq[:eq_ygap], endo[:ygap_t]]  =  1.
       if haskey(m.settings, :regime_eqcond_info)
           if reg >= minimum(keys(get_setting(m, :regime_eqcond_info))) &&
               reg <= maximum(keys(get_setting(m, :regime_eqcond_info))) &&
               haskey(m.settings, :ygap_type)
               if get_setting(m, :ygap_type) in [:smooth_ait_gdp, :smooth_ait_gdp_alt, :flexible_ait, :rw]
                   Thalf  = haskey(get_settings(m), :gdp_Thalf) ? get_setting(m, :gdp_Thalf) : 10.
                   ρ_ygap = exp(log(0.5) / Thalf)

                   Γ0[eq[:eq_ygap], endo[:ygap_t]] = 1.
                   # Γ0[eq[:eq_ygap], endo[:π_t]]    = -1.
                   Γ1[eq[:eq_ygap], endo[:ygap_t]] = ρ_ygap

                   Γ0[eq[:eq_ygap], endo[:y_t]]    = -1.
                   Γ0[eq[:eq_ygap], endo[:z_t]]    = -1.
                   Γ1[eq[:eq_ygap], endo[:y_t]]    = -1.
               end
               if (haskey(get_settings(m), :add_initialize_pgap_ygap_pseudoobs) ?
                   get_setting(m, :add_initialize_pgap_ygap_pseudoobs) : false)
                   Ψ[eq[:eq_ygap], exo[:ygap_sh]] = 1.
               end
           end
       end
       if (haskey(get_settings(m), :add_initialize_pgap_ygap_pseudoobs) ?
                   get_setting(m, :add_initialize_pgap_ygap_pseudoobs) : false)
                   Ψ[eq[:eq_ygap], exo[:ygap_sh]] = 1.
               end
   end

   if haskey(m.settings, :add_rw) ? get_setting(m, :add_rw) : false
       Γ0[eq[:eq_rw], endo[:rw_t]]     =  1.
       Γ0[eq[:eq_Rref], endo[:Rref_t]] =  1.

       if haskey(m.settings, :regime_eqcond_info)
           if reg >= minimum(keys(get_setting(m, :regime_eqcond_info))) &&
               reg <= maximum(keys(get_setting(m, :regime_eqcond_info))) &&
               haskey(m.settings, :Rref_type)

               ρ_rw = haskey(get_settings(m), :ρ_rw) ? get_setting(m, :ρ_rw) : 0.93
               Γ0[eq[:eq_rw], endo[:rw_t]]   = 1.
               Γ1[eq[:eq_rw], endo[:rw_t]]   = ρ_rw

               if get_setting(m, :Rref_type) in [:ait]
                   Thalf  = haskey(get_settings(m), :ait_Thalf) ? get_setting(m, :ait_Thalf) : 10.
                   ρ_pgap = exp(log(0.5) / Thalf)
                   φ      = haskey(get_settings(m), :ait_φ) ? get_setting(m, :ait_φ) : 0.25

                   Γ0[eq[:eq_Rref], endo[:Rref_t]] = 1.
                   Γ0[eq[:eq_Rref], endo[:pgap_t]] = -φ * (1/(1 - ρ_pgap))
                   Γ1[eq[:eq_Rref], endo[:Rref_t]] = 0.
                   C[eq[:eq_Rref]]                 = 0.
                   Γ0[eq[:eq_Rref], endo[:rw_t]]   = -1.

               elseif get_setting(m, :Rref_type) in [:smooth_ait_gdp, :smooth_ait_gdp_alt, :flexible_ait, :rw]

                   ait_Thalf = haskey(get_settings(m), :ait_Thalf) ? get_setting(m, :ait_Thalf) : 10.
                   gdp_Thalf = haskey(get_settings(m), :gdp_Thalf) ? get_setting(m, :gdp_Thalf) : 10.
                   ρ_pgap    = exp(log(0.5) / ait_Thalf)
                   ρ_ygap    = exp(log(0.5) / gdp_Thalf)
                   ρ_smooth  = haskey(get_settings(m), :rw_ρ_smooth) ? get_setting(m, :rw_ρ_smooth) : 0.656
                   φ_π       = haskey(get_settings(m), :rw_φ_π) ? get_setting(m, :rw_φ_π) : 11.13
                   φ_y       = haskey(get_settings(m), :rw_φ_y) ? get_setting(m, :rw_φ_y) : 11.13

                   Γ0[eq[:eq_Rref], endo[:Rref_t]] = 1.
                   Γ1[eq[:eq_Rref], endo[:Rref_t]] = ρ_smooth
                   C[eq[:eq_Rref]]                 = 0.

                   Γ0[eq[:eq_Rref], endo[:pgap_t]]   = -φ_π * (1. - ρ_pgap) * (1. - ρ_smooth) # This is the AIT part
                   Γ0[eq[:eq_Rref], endo[:ygap_t]]   = -φ_y * (1. - ρ_ygap) * (1. - ρ_smooth) # This is the GDP part
               end
           end
       end
   end

   if haskey(m.settings, :add_pgap)
       if get_setting(m, :add_pgap)
           if get_setting(m, :pgap_type) in [:smooth_ait_gdp, :smooth_ait, :ait, :smooth_ait_gdp_alt, :flexible_ait]
               Thalf = haskey(m.settings, :ait_Thalf) ? get_setting(m, :ait_Thalf) : 10
               ρ_pgap = exp(log(0.5) / Thalf)
               Γ0[eq[:eq_pgap], endo[:pgap_t]] = 1.
               Γ0[eq[:eq_pgap], endo[:π_t]]    = -1.
               Γ1[eq[:eq_pgap], endo[:pgap_t]] = ρ_pgap
               if (haskey(get_settings(m), :add_initialize_pgap_ygap_pseudoobs) ?
                   get_setting(m, :add_initialize_pgap_ygap_pseudoobs) : false)
                   Ψ[eq[:eq_pgap], exo[:pgap_sh]] = 1.
               end
           end
       end
   end

   if haskey(m.settings, :add_ygap)
       if get_setting(m, :add_ygap)
           if get_setting(m, :ygap_type) in [:smooth_ait_gdp, :smooth_ait_gdp_alt, :flexible_ait]
               Thalf  = haskey(get_settings(m), :gdp_Thalf) ? get_setting(m, :gdp_Thalf) : 10.
               ρ_ygap = exp(log(0.5) / Thalf)

               Γ0[eq[:eq_ygap], endo[:ygap_t]] = 1.
               # Γ0[eq[:eq_ygap], endo[:π_t]]    = -1.
               Γ1[eq[:eq_ygap], endo[:ygap_t]] = ρ_ygap

               Γ0[eq[:eq_ygap], endo[:y_t]]    = -1.
               Γ0[eq[:eq_ygap], endo[:z_t]]    = -1.
               Γ1[eq[:eq_ygap], endo[:y_t]]    = -1.
               if (haskey(get_settings(m), :add_initialize_pgap_ygap_pseudoobs) ?
                   get_setting(m, :add_initialize_pgap_ygap_pseudoobs) : false)
                   Ψ[eq[:eq_ygap], exo[:ygap_sh]] = 1.
               end
           end
       end
   end

   for para in m.parameters
        if !isempty(para.regimes)
            ModelConstructors.toggle_regime!(para, 1)
        end
    end

    return Γ0, Γ1, C, Ψ, Π
end
