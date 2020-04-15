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
* `C`  (`n_statss54ges` x `1`) is a vector of constants
* `Ψ`  (`n_states` x `n_shocks_exogenous`) holds coefficients of iid shocks.
* `Π`  (`n_states` x `n_states_expectational`) holds coefficients of expectational states.
"""
function eqcond_regimes(m::Model1002)
    Γ0 = Vector{Array}(undef, 2)
    Γ1 = Vector{Array}(undef, 2)
    C = Vector{Array}(undef, 2)
    Ψ = Vector{Array}(undef, 2)
    Π = Vector{Array}(undef, 2)
    for regime=1:2
        endo = m.endogenous_states
        exo  = m.exogenous_shocks
        ex   = m.expected_shocks
        eq   = m.equilibrium_conditions

        Γ0[regime] = zeros(n_states(m), n_states(m))
        Γ1[regime] = zeros(n_states(m), n_states(m))
        C[regime]  = zeros(n_states(m))
        Ψ[regime]  = zeros(n_states(m), n_shocks_exogenous(m))
        Π[regime]  = zeros(n_states(m), n_shocks_expectational(m))

        ### ENDOGENOUS STATES ###

        ### 1. Consumption Euler Equation

        # Sticky prices and wages
        Γ0[regime][eq[:eq_euler], endo[:c_t]]  = 1.
        Γ0[regime][eq[:eq_euler], endo[:R_t]]  = (1 - m[:h]*exp(-m[:z_star]))/(m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))
        Γ0[regime][eq[:eq_euler], endo[:b_t]]  = -1.
        Γ0[regime][eq[:eq_euler], endo[:Eπ_t]] = -(1 - m[:h]*exp(-m[:z_star]))/(m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))
        Γ0[regime][eq[:eq_euler], endo[:z_t]]  = (m[:h]*exp(-m[:z_star]))/(1 + m[:h]*exp(-m[:z_star]))
        Γ0[regime][eq[:eq_euler], endo[:Ec_t]] = -1/(1 + m[:h]*exp(-m[:z_star]))
        Γ0[regime][eq[:eq_euler], endo[:Ez_t]] = -1/(1 + m[:h]*exp(-m[:z_star]))
        Γ0[regime][eq[:eq_euler], endo[:L_t]]  = -(m[:σ_c] - 1)*m[:wl_c]/(m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))
        Γ0[regime][eq[:eq_euler], endo[:EL_t]] = (m[:σ_c] - 1)*m[:wl_c]/(m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))
        Γ1[regime][eq[:eq_euler], endo[:c_t]]  = (m[:h]*exp(-m[:z_star]))/(1 + m[:h]*exp(-m[:z_star]))

        # Flexible prices and wages
        Γ0[regime][eq[:eq_euler_f], endo[:c_f_t]]  = 1.
        Γ0[regime][eq[:eq_euler_f], endo[:r_f_t]]  = (1 - m[:h]*exp(-m[:z_star]))/(m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))
        Γ0[regime][eq[:eq_euler_f], endo[:b_t]]    = -1.
        Γ0[regime][eq[:eq_euler_f], endo[:z_t]]    = (m[:h]*exp(-m[:z_star]))/(1 + m[:h]*exp(-m[:z_star]))
        Γ0[regime][eq[:eq_euler_f], endo[:Ec_f_t]] = -1/(1 + m[:h]*exp(-m[:z_star]))
        Γ0[regime][eq[:eq_euler_f], endo[:Ez_t]]   = -1/(1 + m[:h]*exp(-m[:z_star]))
        Γ0[regime][eq[:eq_euler_f], endo[:L_f_t]]  = -(m[:σ_c] - 1)*m[:wl_c]/(m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))
        Γ0[regime][eq[:eq_euler_f], endo[:EL_f_t]] = (m[:σ_c] - 1)*m[:wl_c]/(m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))
        Γ1[regime][eq[:eq_euler_f], endo[:c_f_t]]  = (m[:h]*exp(-m[:z_star]))/(1 + m[:h]*exp(-m[:z_star]))

        ### 2. Investment Euler Equation

        # Sticky prices and wages
        Γ0[regime][eq[:eq_inv], endo[:qk_t]] = -1/(m[:S′′]*exp(2.0*m[:z_star])*(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star])))
        Γ0[regime][eq[:eq_inv], endo[:i_t]]  = 1.
        Γ0[regime][eq[:eq_inv], endo[:z_t]]  = 1/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        Γ1[regime][eq[:eq_inv], endo[:i_t]]  = 1/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        Γ0[regime][eq[:eq_inv], endo[:Ei_t]] = -m[:β]*exp((1 - m[:σ_c])*m[:z_star])/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        Γ0[regime][eq[:eq_inv], endo[:Ez_t]] = -m[:β]*exp((1 - m[:σ_c])*m[:z_star])/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        Γ0[regime][eq[:eq_inv], endo[:μ_t]]  = -1.

        # Flexible prices and wages
        Γ0[regime][eq[:eq_inv_f], endo[:qk_f_t]] = -1/(m[:S′′]*exp(2*m[:z_star])*(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star])))
        Γ0[regime][eq[:eq_inv_f], endo[:i_f_t]]  = 1.
        Γ0[regime][eq[:eq_inv_f], endo[:z_t]]    = 1/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        Γ1[regime][eq[:eq_inv_f], endo[:i_f_t]]  = 1/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        Γ0[regime][eq[:eq_inv_f], endo[:Ei_f_t]] = -m[:β]*exp((1 - m[:σ_c])*m[:z_star])/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        Γ0[regime][eq[:eq_inv_f], endo[:Ez_t]]   = -m[:β]*exp((1 - m[:σ_c])*m[:z_star])/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        Γ0[regime][eq[:eq_inv_f], endo[:μ_t]]    = -1.

        ### 3. Financial Friction Block

        # Return to capital
        # Sticky prices and wages
        Γ0[regime][eq[:eq_capval], endo[:Rktil_t]] = 1.
        Γ0[regime][eq[:eq_capval], endo[:π_t]]      = -1.
        Γ0[regime][eq[:eq_capval], endo[:rk_t]]     = -m[:r_k_star]/(1 + m[:r_k_star] - m[:δ])
        Γ0[regime][eq[:eq_capval], endo[:qk_t]]     = -(1 - m[:δ])/(1 + m[:r_k_star] - m[:δ])
        Γ1[regime][eq[:eq_capval], endo[:qk_t]]     = -1.

        # Spreads
        # Sticky prices and wages
        Γ0[regime][eq[:eq_spread], endo[:ERtil_k_t]] = 1.
        Γ0[regime][eq[:eq_spread], endo[:R_t]]       = -1.
        Γ0[regime][eq[:eq_spread], endo[:b_t]]       = (m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))/(1 - m[:h]*exp(-m[:z_star]))
        Γ0[regime][eq[:eq_spread], endo[:qk_t]]      = -m[:ζ_spb]
        Γ0[regime][eq[:eq_spread], endo[:kbar_t]]    = -m[:ζ_spb]
        Γ0[regime][eq[:eq_spread], endo[:n_t]]       = m[:ζ_spb]
        Γ0[regime][eq[:eq_spread], endo[:σ_ω_t]]     = -1.
        Γ0[regime][eq[:eq_spread], endo[:μ_e_t]]     = -1.

        # Flexible prices and wages
        Γ0[regime][eq[:eq_spread_f], endo[:ERktil_f_t]] = 1.
        Γ0[regime][eq[:eq_spread_f], endo[:r_f_t]]       = -1.
        Γ0[regime][eq[:eq_spread_f], endo[:b_t]]       = (m[:σ_c]*(1 + m[:h]*exp(-m[:z_star])))/(1 - m[:h]*exp(-m[:z_star]))
        Γ0[regime][eq[:eq_spread_f], endo[:qk_f_t]]      = -m[:ζ_spb]
        Γ0[regime][eq[:eq_spread_f], endo[:kbar_f_t]]    = -m[:ζ_spb]
        Γ0[regime][eq[:eq_spread_f], endo[:n_f_t]]       = m[:ζ_spb]
        Γ0[regime][eq[:eq_spread_f], endo[:σ_ω_t]]     = -1.
        Γ0[regime][eq[:eq_spread_f], endo[:μ_e_t]]     = -1.

        # n evol
        # Sticky prices and wages
        Γ0[regime][eq[:eq_nevol], endo[:n_t]]      = 1.
        Γ0[regime][eq[:eq_nevol], endo[:γ_t]]      = -1.
        Γ0[regime][eq[:eq_nevol], endo[:z_t]]      = m[:γ_star]*m[:vstar]/m[:nstar]
        Γ0[regime][eq[:eq_nevol], endo[:Rktil_t]] = -m[:ζ_nRk]
        Γ0[regime][eq[:eq_nevol], endo[:π_t]]      = (m[:ζ_nRk] - m[:ζ_nR])
        Γ1[regime][eq[:eq_nevol], endo[:σ_ω_t]]    = -m[:ζ_nσ_ω]/m[:ζ_spσ_ω]
        Γ1[regime][eq[:eq_nevol], endo[:μ_e_t]]    = -m[:ζ_nμ_e]/m[:ζ_spμ_e]
        Γ1[regime][eq[:eq_nevol], endo[:qk_t]]     = m[:ζ_nqk]
        Γ1[regime][eq[:eq_nevol], endo[:kbar_t]]   = m[:ζ_nqk]
        Γ1[regime][eq[:eq_nevol], endo[:n_t]]      = m[:ζ_nn]
        Γ1[regime][eq[:eq_nevol], endo[:R_t]]      = -m[:ζ_nR]
        Γ1[regime][eq[:eq_nevol], endo[:b_t]]      = m[:ζ_nR]*((m[:σ_c]*(1.0+m[:h]*exp(-m[:z_star])))/(1.0-m[:h]*exp(-m[:z_star])))

        # Flexible prices and wages
        Γ0[regime][eq[:eq_nevol_f], endo[:n_f_t]]      = 1.
        Γ0[regime][eq[:eq_nevol_f], endo[:z_t]]      = m[:γ_star]*m[:vstar]/m[:nstar]
        Γ0[regime][eq[:eq_nevol_f], endo[:rktil_f_t]] = -m[:ζ_nRk]
        Γ1[regime][eq[:eq_nevol_f], endo[:σ_ω_t]]    = -m[:ζ_nσ_ω]/m[:ζ_spσ_ω]
        Γ1[regime][eq[:eq_nevol_f], endo[:μ_e_t]]    = -m[:ζ_nμ_e]/m[:ζ_spμ_e]
        Γ1[regime][eq[:eq_nevol_f], endo[:qk_f_t]]     = m[:ζ_nqk]
        Γ1[regime][eq[:eq_nevol_f], endo[:kbar_f_t]]   = m[:ζ_nqk]
        Γ1[regime][eq[:eq_nevol_f], endo[:n_f_t]]      = m[:ζ_nn]
        Γ1[regime][eq[:eq_nevol_f], endo[:r_f_t]]      = -m[:ζ_nR]
        Γ1[regime][eq[:eq_nevol_f], endo[:b_t]]      = m[:ζ_nR]*((m[:σ_c]*(1.0+m[:h]*exp(-m[:z_star])))/(1.0-m[:h]*exp(-m[:z_star])))

        # Flexible prices and wages - ASSUME NO FINANCIAL FRICTIONS
        Γ0[regime][eq[:eq_capval_f], endo[:rktil_f_t]] = 1.
        Γ0[regime][eq[:eq_capval_f], endo[:rk_f_t]]     = -m[:r_k_star]/(m[:r_k_star]+1-m[:δ])
        Γ0[regime][eq[:eq_capval_f], endo[:qk_f_t]]     = -(1-m[:δ])/(m[:r_k_star]+1-m[:δ])
        Γ1[regime][eq[:eq_capval_f], endo[:qk_f_t]]     = -1.

        ### 4. Aggregate Production Function

        # Sticky prices and wages
        Γ0[regime][eq[:eq_output], endo[:y_t]] =  1.
        if subspec(m) in ["ss47", "ss48"] && regime == 2
            Γ0[regime][eq[:eq_output], endo[:k_t]] = -m[:Φ_r2]*m[:α]
            Γ0[regime][eq[:eq_output], endo[:L_t]] = -m[:Φ_r2]*(1 - m[:α])
        else
            Γ0[regime][eq[:eq_output], endo[:k_t]] = -m[:Φ]*m[:α]
            Γ0[regime][eq[:eq_output], endo[:L_t]] = -m[:Φ]*(1 - m[:α])
        end

        # Flexible prices and wages
        Γ0[regime][eq[:eq_output_f], endo[:y_f_t]] =  1.
        if subspec(m) in ["ss47", "ss48"] && regime == 2
            Γ0[regime][eq[:eq_output_f], endo[:k_f_t]] = -m[:Φ_r2]*m[:α]
            Γ0[regime][eq[:eq_output_f], endo[:L_f_t]] = -m[:Φ_r2]*(1 - m[:α])
        else
            Γ0[regime][eq[:eq_output_f], endo[:k_f_t]] = -m[:Φ]*m[:α]
            Γ0[regime][eq[:eq_output_f], endo[:L_f_t]] = -m[:Φ]*(1 - m[:α])
        end

        ### 5. Capital Utilization

        # Sticky prices and wages
        Γ0[regime][eq[:eq_caputl], endo[:k_t]]    =  1.
        Γ1[regime][eq[:eq_caputl], endo[:kbar_t]] =  1.
        Γ0[regime][eq[:eq_caputl], endo[:z_t]]    = 1.
        Γ0[regime][eq[:eq_caputl], endo[:u_t]]    = -1.

        # Flexible prices and wages
        Γ0[regime][eq[:eq_caputl_f], endo[:k_f_t]]    =  1.
        Γ1[regime][eq[:eq_caputl_f], endo[:kbar_f_t]] =  1.
        Γ0[regime][eq[:eq_caputl_f], endo[:z_t]]      = 1.
        Γ0[regime][eq[:eq_caputl_f], endo[:u_f_t]]    = -1.

        ### 6. Rental Rate of Capital

        # Sticky prices and wages
        Γ0[regime][eq[:eq_capsrv], endo[:u_t]]  = 1.
        Γ0[regime][eq[:eq_capsrv], endo[:rk_t]] = -(1 - m[:ppsi])/m[:ppsi]

        # Flexible prices and wages
        Γ0[regime][eq[:eq_capsrv_f], endo[:u_f_t]]  = 1.
        Γ0[regime][eq[:eq_capsrv_f], endo[:rk_f_t]] = -(1 - m[:ppsi])/m[:ppsi]

        ### 7. Evolution of Capital

        # Sticky prices and wages
        Γ0[regime][eq[:eq_capev], endo[:kbar_t]] = 1.
        Γ1[regime][eq[:eq_capev], endo[:kbar_t]] = 1 - m[:istar]/m[:kbarstar]
        Γ0[regime][eq[:eq_capev], endo[:z_t]]    = 1 - m[:istar]/m[:kbarstar]
        Γ0[regime][eq[:eq_capev], endo[:i_t]]    = -m[:istar]/m[:kbarstar]
        Γ0[regime][eq[:eq_capev], endo[:μ_t]]    = -m[:istar]*m[:S′′]*exp(2*m[:z_star])*(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))/m[:kbarstar]

        # Flexible prices and wages
        Γ0[regime][eq[:eq_capev_f], endo[:kbar_f_t]] = 1.
        Γ1[regime][eq[:eq_capev_f], endo[:kbar_f_t]] = 1 - m[:istar]/m[:kbarstar]
        Γ0[regime][eq[:eq_capev_f], endo[:z_t]]      = 1 - m[:istar]/m[:kbarstar]
        Γ0[regime][eq[:eq_capev_f], endo[:i_f_t]]    = -m[:istar]/m[:kbarstar]
        Γ0[regime][eq[:eq_capev_f], endo[:μ_t]]      = -m[:istar]*m[:S′′]*exp(2*m[:z_star])*(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))/m[:kbarstar]

        ### 8. Price Markup

        # Sticky prices and wages
        Γ0[regime][eq[:eq_mkupp], endo[:mc_t]] =  1.
        Γ0[regime][eq[:eq_mkupp], endo[:w_t]]  = -1.
        Γ0[regime][eq[:eq_mkupp], endo[:L_t]]  = -m[:α]
        Γ0[regime][eq[:eq_mkupp], endo[:k_t]]  =  m[:α]

        # Flexible prices and wages
        Γ0[regime][eq[:eq_mkupp_f], endo[:w_f_t]] = 1.
        Γ0[regime][eq[:eq_mkupp_f], endo[:L_f_t]] =  m[:α]
        Γ0[regime][eq[:eq_mkupp_f], endo[:k_f_t]] =  -m[:α]

        ### 9. Phillips Curve

        # Sticky prices and wages
        Γ0[regime][eq[:eq_phlps], endo[:π_t]]  = 1.
#        @show subspec(m), regime
        if subspec(m) in ["ss21", "ss22", "ss25", "ss26", "ss28", "ss29", "ss41", "ss42"] && regime == 2
            Γ0[regime][eq[:eq_phlps], endo[:mc_t]] =  -((1 - m[:ζ_p_r2]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))*
                                                        (1 - m[:ζ_p_r2]))/(m[:ζ_p_r2]*((m[:Φ]- 1)*m[:ϵ_p] + 1))/(1 + m[:ι_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        elseif subspec(m) in ["ss45", "ss46"] && regime == 2
            Γ0[regime][eq[:eq_phlps], endo[:mc_t]] =  -((1 - m[:ζ_p_r2]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))*
                                                        (1 - m[:ζ_p_r2]))/(m[:ζ_p_r2]*((m[:Φ]- 1)*m[:ϵ_p] + 1))/(1 + m[:ι_p_r2]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        elseif subspec(m) in ["ss47", "ss48"] && regime == 2
            Γ0[regime][eq[:eq_phlps], endo[:mc_t]] =  -((1 - m[:ζ_p_r2]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))*
                                                        (1 - m[:ζ_p_r2]))/(m[:ζ_p_r2]*((m[:Φ_r2]- 1)*m[:ϵ_p] + 1))/(1 + m[:ι_p_r2]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        else
            Γ0[regime][eq[:eq_phlps], endo[:mc_t]] =  -((1 - m[:ζ_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))*
                                                        (1 - m[:ζ_p]))/(m[:ζ_p]*((m[:Φ]- 1)*m[:ϵ_p] + 1))/(1 + m[:ι_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        end

        if subspec(m) in ["ss45", "ss46", "ss47", "ss48"] && regime == 2
            Γ1[regime][eq[:eq_phlps], endo[:π_t]]  = m[:ι_p_r2]/(1 + m[:ι_p_r2]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
            Γ0[regime][eq[:eq_phlps], endo[:Eπ_t]] = -m[:β]*exp((1 - m[:σ_c])*m[:z_star])/(1 + m[:ι_p_r2]*m[:β]*
                                                                                           exp((1 - m[:σ_c])*m[:z_star]))
        else
            Γ1[regime][eq[:eq_phlps], endo[:π_t]]  = m[:ι_p]/(1 + m[:ι_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
            Γ0[regime][eq[:eq_phlps], endo[:Eπ_t]] = -m[:β]*exp((1 - m[:σ_c])*m[:z_star])/(1 + m[:ι_p]*m[:β]*
                                                                                           exp((1 - m[:σ_c])*m[:z_star]))
        end

        # Comment out for counterfactual with no price mark up shock
        if subspec(m) == "ss20"
            # Re-scale markup shock so that it is the shock corresponds to the size of the real
            # markup shock. We use this spec for counterfactual analysis when using different
            # parameter values for parameters governing κ, e.g. ζ_p.
            #
            # The standard implementation scales the actual markup shock, say
            # markup_t, by κ, i.e. λ_f = κ * markup_t.
            # The real level and standard deviation of markup_t are λ_f / kappa and
            # σ_λ_f / κ, respectively. The other parameters are accurately estimated
            # in the ARMA(1,1) governing λ_f.
            #
            # To rescale, we set the coefficient on λ_f to be κnum / κden,
            # where κnum = κ (i.e. it is the correct kappa), but
            # κden fixes the value of specific parameters and retrieves those desired values
            # from the settings, which are provided through custom_settings during
            # instantiation of a model object. In this way,
            # we can back out the correct markup shock after estimating
            # on the more numerically stable implementation we usually use.
            κnum = ((1 - m[:ζ_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))*
                    (1 - m[:ζ_p]))/(m[:ζ_p]*((m[:Φ]- 1)*m[:ϵ_p] + 1))/(1 + m[:ι_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))         # kappa numerator
            fix_ζ_p = get_setting(m, :fix_ζ_p)
            κden = ((1 - fix_ζ_p*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))*
                    (1 - fix_ζ_p))/(fix_ζ_p*((m[:Φ]- 1)*m[:ϵ_p] + 1))/(1 + m[:ι_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))         # kappa denominator

            Γ0[regime][eq[:eq_phlps], endo[:λ_f_t]] = -κnum / κden
        #=elseif subspec(m) in ["ss21", "ss22", "ss25", "ss26", "ss28", "ss29", "ss41", "ss42"]
            if regime == 1
                Γ0[regime][eq[:eq_phlps], endo[:λ_f_t]] = -1.
            elseif regime == 2
                κnum = ((1 - m[:ζ_p_r2]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))*
                        (1 - m[:ζ_p_r2]))/(m[:ζ_p_r2]*((m[:Φ]- 1)*m[:ϵ_p] + 1))/(1 + m[:ι_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))         # kappa numerator
                fix_ζ_p = m[:ζ_p]
                κden = ((1 - fix_ζ_p*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))*
                        (1 - fix_ζ_p))/(fix_ζ_p*((m[:Φ]- 1)*m[:ϵ_p] + 1))/(1 + m[:ι_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))         # kappa denominator
                Γ0[regime][eq[:eq_phlps], endo[:λ_f_t]] = -κnum / κden
            end=#
        else
            Γ0[regime][eq[:eq_phlps], endo[:λ_f_t]] = -1.
        end

        # Flexible prices and wages not necessary

        ### 10. Rental Rate of Capital

        # Sticky prices and wages
        Γ0[regime][eq[:eq_caprnt], endo[:rk_t]] = 1.
        Γ0[regime][eq[:eq_caprnt], endo[:k_t]]  = 1.
        Γ0[regime][eq[:eq_caprnt], endo[:L_t]]  = -1.
        Γ0[regime][eq[:eq_caprnt], endo[:w_t]]  = -1.

        # Flexible prices and wages
        Γ0[regime][eq[:eq_caprnt_f], endo[:rk_f_t]] = 1.
        Γ0[regime][eq[:eq_caprnt_f], endo[:k_f_t]]  = 1.
        Γ0[regime][eq[:eq_caprnt_f], endo[:L_f_t]]  = -1.
        Γ0[regime][eq[:eq_caprnt_f], endo[:w_f_t]]  = -1.

        ### 11. Marginal Substitution

        # Sticky prices and wages
        Γ0[regime][eq[:eq_msub], endo[:μ_ω_t]] = 1.
        Γ0[regime][eq[:eq_msub], endo[:L_t]]   = m[:ν_l]
        Γ0[regime][eq[:eq_msub], endo[:c_t]]   = 1/(1 - m[:h]*exp(-m[:z_star]))
        Γ1[regime][eq[:eq_msub], endo[:c_t]]   = m[:h]*exp(-m[:z_star])/(1 - m[:h]*exp(-m[:z_star]))
        Γ0[regime][eq[:eq_msub], endo[:z_t]]   = m[:h]*exp(-m[:z_star]) /(1 - m[:h]*exp(-m[:z_star]))
        Γ0[regime][eq[:eq_msub], endo[:w_t]]   = -1.

        # Flexible prices and wages
        Γ0[regime][eq[:eq_msub_f], endo[:w_f_t]] = -1.
        Γ0[regime][eq[:eq_msub_f], endo[:L_f_t]] = m[:ν_l]
        Γ0[regime][eq[:eq_msub_f], endo[:c_f_t]] = 1/(1 - m[:h]*exp(-m[:z_star]))
        Γ1[regime][eq[:eq_msub_f], endo[:c_f_t]] = m[:h]*exp(-m[:z_star])/(1 - m[:h]*exp(-m[:z_star]))
        Γ0[regime][eq[:eq_msub_f], endo[:z_t]]   = m[:h]*exp(-m[:z_star])/(1 - m[:h]*exp(-m[:z_star]))

        ### 12. Evolution of Wages

        # Sticky prices and wages
        Γ0[regime][eq[:eq_wage], endo[:w_t]]   = 1
        Γ0[regime][eq[:eq_wage], endo[:μ_ω_t]] = (1 - m[:ζ_w]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))*
(1 - m[:ζ_w])/(m[:ζ_w]*((m[:λ_w] - 1)*m[:ϵ_w] + 1))/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        Γ0[regime][eq[:eq_wage], endo[:π_t]]   = (1 + m[:ι_w]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        Γ1[regime][eq[:eq_wage], endo[:w_t]]   = 1/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        Γ0[regime][eq[:eq_wage], endo[:z_t]]   = 1/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        Γ1[regime][eq[:eq_wage], endo[:π_t]]   = m[:ι_w]/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        Γ0[regime][eq[:eq_wage], endo[:Ew_t]]  = -m[:β]*exp((1 - m[:σ_c])*m[:z_star])/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        Γ0[regime][eq[:eq_wage], endo[:Ez_t]]  = -m[:β]*exp((1 - m[:σ_c])*m[:z_star])/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        Γ0[regime][eq[:eq_wage], endo[:Eπ_t]]  = -m[:β]*exp((1 - m[:σ_c])*m[:z_star])/(1 + m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        Γ0[regime][eq[:eq_wage], endo[:λ_w_t]] = -1.

        # Flexible prices and wages not necessary

        ### 13. Monetary Policy Rule

        # Sticky prices and wages
        Γ0[regime][eq[:eq_mp], endo[:R_t]]      = 1.
        if subspec(m) in ["ss23", "ss24", "ss25", "ss26", "ss28", "ss29", "ss43", "ss44"] && regime == 2
            Γ1[regime][eq[:eq_mp], endo[:R_t]]      = m[:ρ_r2]
            Γ0[regime][eq[:eq_mp], endo[:π_t]]      = -(1 - m[:ρ_r2])*m[:ψ1_r2]
            Γ0[regime][eq[:eq_mp], endo[:π_star_t]] = (1 - m[:ρ_r2])*m[:ψ1_r2]
            Γ0[regime][eq[:eq_mp], endo[:y_t]]      = -(1 - m[:ρ_r2])*m[:ψ2_r2] - m[:ψ3_r2 ]
            Γ0[regime][eq[:eq_mp], endo[:y_f_t]]    = (1 - m[:ρ_r2])*m[:ψ2_r2] + m[:ψ3_r2]
            Γ1[regime][eq[:eq_mp], endo[:y_t]]      = -m[:ψ3_r2]
            Γ1[regime][eq[:eq_mp], endo[:y_f_t]]    = m[:ψ3_r2]
        else
            Γ1[regime][eq[:eq_mp], endo[:R_t]]      = m[:ρ]
            Γ0[regime][eq[:eq_mp], endo[:π_t]]      = -(1 - m[:ρ])*m[:ψ1]
            Γ0[regime][eq[:eq_mp], endo[:π_star_t]] = (1 - m[:ρ])*m[:ψ1]
            Γ0[regime][eq[:eq_mp], endo[:y_t]]      = -(1 - m[:ρ])*m[:ψ2] - m[:ψ3]
            Γ0[regime][eq[:eq_mp], endo[:y_f_t]]    = (1 - m[:ρ])*m[:ψ2] + m[:ψ3]
            Γ1[regime][eq[:eq_mp], endo[:y_t]]      = -m[:ψ3]
            Γ1[regime][eq[:eq_mp], endo[:y_f_t]]    = m[:ψ3]
        end
        Γ0[regime][eq[:eq_mp], endo[:rm_t]]     = -1.

        # Flexible prices and wages not necessary

        ### 14. Resource Constraint

        # Sticky prices and wages
        Γ0[regime][eq[:eq_res], endo[:y_t]] = 1.
        Γ0[regime][eq[:eq_res], endo[:g_t]] = -m[:g_star]
        Γ0[regime][eq[:eq_res], endo[:c_t]] = -m[:cstar]/m[:ystar]
        Γ0[regime][eq[:eq_res], endo[:i_t]] = -m[:istar]/m[:ystar]
        Γ0[regime][eq[:eq_res], endo[:u_t]] = -m[:r_k_star]*m[:kstar]/m[:ystar]

        # Flexible prices and wages
        Γ0[regime][eq[:eq_res_f], endo[:y_f_t]] = 1.
        Γ0[regime][eq[:eq_res_f], endo[:g_t]]   = -m[:g_star]
        Γ0[regime][eq[:eq_res_f], endo[:c_f_t]] = -m[:cstar]/m[:ystar]
        Γ0[regime][eq[:eq_res_f], endo[:i_f_t]] = -m[:istar]/m[:ystar]
        Γ0[regime][eq[:eq_res_f], endo[:u_f_t]] = -m[:r_k_star]*m[:kstar]/m[:ystar]

        ### 15. Extra States
        # These aren't strictly necessary, but they track lags or simplify the equations

        # π_t1
        Γ0[regime][eq[:eq_π1], endo[:π_t1]] = 1.
        Γ1[regime][eq[:eq_π1], endo[:π_t]]  = 1.

        # π_t2
        Γ0[regime][eq[:eq_π2], endo[:π_t2]] = 1.
        Γ1[regime][eq[:eq_π2], endo[:π_t1]] = 1.

        # π_a
        Γ0[regime][eq[:eq_π_a], endo[:π_a_t]] = 1.
        Γ0[regime][eq[:eq_π_a], endo[:π_t]]   = -1.
        Γ0[regime][eq[:eq_π_a], endo[:π_t1]]  = -1.
        Γ0[regime][eq[:eq_π_a], endo[:π_t2]]  = -1.
        Γ1[regime][eq[:eq_π_a], endo[:π_t2]]  = 1.

        # Rt1
        Γ0[regime][eq[:eq_Rt1], endo[:R_t1]] = 1.
        Γ1[regime][eq[:eq_Rt1], endo[:R_t]]  = 1.

        # Ez_t
        Γ0[regime][eq[:eq_Ez], endo[:Ez_t]]   = 1.
        Γ0[regime][eq[:eq_Ez], endo[:ztil_t]]   = -(m[:ρ_ztil]-1)/(1-m[:α])
        Γ0[regime][eq[:eq_Ez], endo[:zp_t]]   = -m[:ρ_z_p]

        ### EXOGENOUS SHOCKS ###

        # Neutral technology
        Γ0[regime][eq[:eq_z], endo[:z_t]]    = 1.
        Γ0[regime][eq[:eq_z], endo[:ztil_t]]  = -1 / (1 - m[:α])
        Γ1[regime][eq[:eq_z], endo[:ztil_t]]  = -1 / (1 - m[:α])
        Γ0[regime][eq[:eq_z], endo[:zp_t]]   = -1.

        Γ0[regime][eq[:eq_ztil], endo[:ztil_t]] = 1.
        Γ1[regime][eq[:eq_ztil], endo[:ztil_t]] = m[:ρ_ztil]
        Ψ[regime][eq[:eq_ztil], exo[:ztil_sh]]     = 1.

        if subspec(m) == "ss60"
            # Ez_t
            Γ0[regime][eq[:eq_Ez], endo[:ziid_t]]   = -(m[:ρ_ziid]-1)/(1-m[:α])

            # Neutral technology
            Γ0[regime][eq[:eq_z], endo[:ziid_t]]  = -1 / (1 - m[:α])
            Γ1[regime][eq[:eq_z], endo[:ziid_t]]  = -1 / (1 - m[:α])

            # AR(1) for ziid
            Γ0[regime][eq[:eq_ziid], endo[:ziid_t]] = 1.
            Γ1[regime][eq[:eq_ziid], endo[:ziid_t]] = m[:ρ_ziid]
            Ψ[regime][eq[:eq_ziid], exo[:ziid_sh]]     = 1.
        end

        # Long-run changes to productivity
        Γ0[regime][eq[:eq_zp], endo[:zp_t]] = 1.
        Γ1[regime][eq[:eq_zp], endo[:zp_t]] = m[:ρ_z_p]
        Ψ[regime][eq[:eq_zp], exo[:zp_sh]]  = 1.

        # Government spending
        Γ0[regime][eq[:eq_g], endo[:g_t]] = 1.
        Γ1[regime][eq[:eq_g], endo[:g_t]] = m[:ρ_g]
        Ψ[regime][eq[:eq_g], exo[:g_sh]]  = 1.
        Ψ[regime][eq[:eq_g], exo[:ztil_sh]]  = m[:η_gz]

        # Asset shock
        Γ0[regime][eq[:eq_b], endo[:b_t]] = 1.
        Γ1[regime][eq[:eq_b], endo[:b_t]] = m[:ρ_b]
        Ψ[regime][eq[:eq_b], exo[:b_sh]]  = 1.

        if subspec(m) == "ss60"
            # iid shock
            Γ0[regime][eq[:eq_biid], endo[:biid_t]] = 1.
            Γ1[regime][eq[:eq_biid], endo[:biid_t]] = m[:ρ_biid]
            Ψ[regime][eq[:eq_biid], exo[:biid_sh]]  = 1.

            # Add to AR(1) shock
            Γ0[regime][eq[:eq_b], endo[:biid_t]] = -1.
            Γ1[regime][eq[:eq_b], endo[:biid_t]] = -1.
        end

        # Investment-specific technology
        Γ0[regime][eq[:eq_μ], endo[:μ_t]] = 1.
        Γ1[regime][eq[:eq_μ], endo[:μ_t]] = m[:ρ_μ]
        Ψ[regime][eq[:eq_μ], exo[:μ_sh]]  = 1.

        # Price mark-up shock
        Γ0[regime][eq[:eq_λ_f], endo[:λ_f_t]]  = 1.
        Γ1[regime][eq[:eq_λ_f], endo[:λ_f_t]]  = m[:ρ_λ_f]
        Γ1[regime][eq[:eq_λ_f], endo[:λ_f_t1]] = -m[:η_λ_f]
        Ψ[regime][eq[:eq_λ_f], exo[:λ_f_sh]]   = 1.

        Γ0[regime][eq[:eq_λ_f1], endo[:λ_f_t1]] = 1.
        Ψ[regime][eq[:eq_λ_f1], exo[:λ_f_sh]]   = 1.

        # Wage mark-up shock
        Γ0[regime][eq[:eq_λ_w], endo[:λ_w_t]]  = 1.
        Γ1[regime][eq[:eq_λ_w], endo[:λ_w_t]]  = m[:ρ_λ_w]
        Γ1[regime][eq[:eq_λ_w], endo[:λ_w_t1]] = -m[:η_λ_w]
        Ψ[regime][eq[:eq_λ_w], exo[:λ_w_sh]]   = 1.

        Γ0[regime][eq[:eq_λ_w1], endo[:λ_w_t1]] = 1.
        Ψ[regime][eq[:eq_λ_w1], exo[:λ_w_sh]]   = 1.

        # Monetary policy shock
        Γ0[regime][eq[:eq_rm], endo[:rm_t]] = 1.
        Γ1[regime][eq[:eq_rm], endo[:rm_t]] = m[:ρ_rm]
        Ψ[regime][eq[:eq_rm], exo[:rm_sh]]  = 1.

        ### Financial frictions

        # Standard deviation of capital shock to entrepreneurs
        Γ0[regime][eq[:eq_σ_ω], endo[:σ_ω_t]] = 1.
        Γ1[regime][eq[:eq_σ_ω], endo[:σ_ω_t]] = m[:ρ_σ_w]
        Ψ[regime][eq[:eq_σ_ω], exo[:σ_ω_sh]]  = 1.

        if subspec(m) == "ss60"
            # iid shock
            Γ0[regime][eq[:eq_σ_ωiid], endo[:σ_ωiid_t]] = 1.
            Γ1[regime][eq[:eq_σ_ωiid], endo[:σ_ωiid_t]] = m[:ρ_σ_ωiid]
            Ψ[regime][eq[:eq_σ_ωiid], exo[:σ_ωiid_sh]]  = 1.

            # Add to AR(1) shock
            Γ0[regime][eq[:eq_σ_ω], endo[:σ_ωiid_t]] = -1.
            Γ1[regime][eq[:eq_σ_ω], endo[:σ_ωiid_t]] = -1.
        end

        # Exogenous bankruptcy costs
        Γ0[regime][eq[:eq_μ_e], endo[:μ_e_t]] = 1.
        Γ1[regime][eq[:eq_μ_e], endo[:μ_e_t]] = m[:ρ_μ_e]
        Ψ[regime][eq[:eq_μ_e], exo[:μ_e_sh]]  = 1.

        # Fraction of entrepreneurs surviving period t
        Γ0[regime][eq[:eq_γ], endo[:γ_t]] = 1.
        Γ1[regime][eq[:eq_γ], endo[:γ_t]] = m[:ρ_γ]
        Ψ[regime][eq[:eq_γ], exo[:γ_sh]]  = 1.

        # Long-term inflation expectations
        Γ0[regime][eq[:eq_π_star], endo[:π_star_t]] = 1.
        Γ1[regime][eq[:eq_π_star], endo[:π_star_t]] = m[:ρ_π_star]
        Ψ[regime][eq[:eq_π_star], exo[:π_star_sh]]  = 1.

        # Anticipated policy shocks
        if n_mon_anticipated_shocks(m) > 0

            # This section adds the anticipated shocks. There is one state for all the
            # anticipated shocks that will hit in a given period (i.e. rm_tl2 holds those that
            # will hit in two periods), and the equations are set up so that rm_tl2 last period
            # will feed into rm_tl1 this period (and so on for other numbers), and last period's
            # rm_tl1 will feed into the rm_t process (and affect the Taylor Rule this period).

            Γ1[regime][eq[:eq_rm], endo[:rm_tl1]]   = 1.
            Γ0[regime][eq[:eq_rml1], endo[:rm_tl1]] = 1.
            Ψ[regime][eq[:eq_rml1], exo[:rm_shl1]]  = 1.

            if n_mon_anticipated_shocks(m) > 1
                for i = 2:n_mon_anticipated_shocks(m)
                    Γ1[regime][eq[Symbol("eq_rml$(i-1)")], endo[Symbol("rm_tl$i")]] = 1.
                    Γ0[regime][eq[Symbol("eq_rml$i")], endo[Symbol("rm_tl$i")]]     = 1.
                    Ψ[regime][eq[Symbol("eq_rml$i")], exo[Symbol("rm_shl$i")]]      = 1.
                end
            end
        end

        for (key, val) in get_setting(m, :antshocks)
            an_eq_mapping = get_setting(m, :ant_eq_mapping)
            if val > 0
            # This section adds the anticipated shocks. There is one state for all the
                # anticipated shocks that will hit in a given period (i.e. rm_tl2 holds those that
                # will hit in two periods), and the equations are set up so that rm_tl2 last period
                # will feed into rm_tl1 this period (and so on for other numbers), and last period's
                # rm_tl1 will feed into the rm_t process (and affect the Taylor Rule this period).

                Γ1[regime][eq[:eq_ztil], endo[:z_tl1]]   = 1.
                Γ0[regime][eq[:eq_zl1], endo[:z_tl1]] = 1.
                Ψ[regime][eq[:eq_zl1], exo[:z_shl1]]  = 1.

                # Ez_t
                Γ0[regime][eq[:eq_Ez], endo[:z_tl1]]  = -1 / (1 - m[:α]) # note z_tl1 is a sum of all shocks that will hit next period.
                                                                         # so this is the only required line

                # Same thing as above for z_p is required, and more generally for any Ez equations w/anticipated shocks
                # Γ0[regime][eq[:eq_Ez], endo[:zp_tl1]]   = -1
                if val > 1
                for i = 2:val
                    Γ1[regime][eq[Symbol("eq_zl$(i-1)")], endo[Symbol("z_tl$i")]] = 1.
                    Γ0[regime][eq[Symbol("eq_zl$i")], endo[Symbol("z_tl$i")]]     = 1.
                    Ψ[regime][eq[Symbol("eq_zl$i")], exo[Symbol("z_shl$i")]]      = 1.
                end
                end
            end
        end

        ### EXPECTATION ERRORS ###

        ### E(c)

        # Sticky prices and wages
        Γ0[regime][eq[:eq_Ec], endo[:c_t]]  = 1.
        Γ1[regime][eq[:eq_Ec], endo[:Ec_t]] = 1.
        Π[regime][eq[:eq_Ec], ex[:Ec_sh]]   = 1.

        # Flexible prices and wages
        Γ0[regime][eq[:eq_Ec_f], endo[:c_f_t]]  = 1.
        Γ1[regime][eq[:eq_Ec_f], endo[:Ec_f_t]] = 1.
        Π[regime][eq[:eq_Ec_f], ex[:Ec_f_sh]]   = 1.

        ### E(q)

        # Sticky prices and wages
        Γ0[regime][eq[:eq_Eqk], endo[:qk_t]]  = 1.
        Γ1[regime][eq[:eq_Eqk], endo[:Eqk_t]] = 1.
        Π[regime][eq[:eq_Eqk], ex[:Eqk_sh]]   = 1.

        # Flexible prices and wages
        Γ0[regime][eq[:eq_Eqk_f], endo[:qk_f_t]]  = 1.
        Γ1[regime][eq[:eq_Eqk_f], endo[:Eqk_f_t]] = 1.
        Π[regime][eq[:eq_Eqk_f], ex[:Eqk_f_sh]]   = 1.

        ### E(i)

        # Sticky prices and wages
        Γ0[regime][eq[:eq_Ei], endo[:i_t]]  = 1.
        Γ1[regime][eq[:eq_Ei], endo[:Ei_t]] = 1.
        Π[regime][eq[:eq_Ei], ex[:Ei_sh]]   = 1.

        # Flexible prices and wages
        Γ0[regime][eq[:eq_Ei_f], endo[:i_f_t]]  = 1.
        Γ1[regime][eq[:eq_Ei_f], endo[:Ei_f_t]] = 1.
        Π[regime][eq[:eq_Ei_f], ex[:Ei_f_sh]]   = 1.

        ### E(π)

        # Sticky prices and wages
        Γ0[regime][eq[:eq_Eπ], endo[:π_t]]  = 1.
        Γ1[regime][eq[:eq_Eπ], endo[:Eπ_t]] = 1.
        Π[regime][eq[:eq_Eπ], ex[:Eπ_sh]]   = 1.

        ### E(l)

        # Sticky prices and wages
        Γ0[regime][eq[:eq_EL], endo[:L_t]]  = 1.
        Γ1[regime][eq[:eq_EL], endo[:EL_t]] = 1.
        Π[regime][eq[:eq_EL], ex[:EL_sh]]   = 1.

        # Flexible prices and wages
        Γ0[regime][eq[:eq_EL_f], endo[:L_f_t]]  = 1.
        Γ1[regime][eq[:eq_EL_f], endo[:EL_f_t]] = 1.
        Π[regime][eq[:eq_EL_f], ex[:EL_f_sh]]   = 1.

        ## # E(rk)

        # Sticky prices and wages
        Γ0[regime][eq[:eq_Erk], endo[:rk_t]]  = 1.
        Γ1[regime][eq[:eq_Erk], endo[:Erk_t]] = 1.
        Π[regime][eq[:eq_Erk], ex[:Erk_sh]]   = 1.

        # Flexible prices and wages
        Γ0[regime][eq[:eq_Erktil_f], endo[:rktil_f_t]]  = 1.
        Γ1[regime][eq[:eq_Erktil_f], endo[:ERktil_f_t]] = 1.
        Π[regime][eq[:eq_Erktil_f], ex[:Erktil_f_sh]]   = 1.

        ### E(w)

        # Sticky prices and wages
        Γ0[regime][eq[:eq_Ew], endo[:w_t]]  = 1.
        Γ1[regime][eq[:eq_Ew], endo[:Ew_t]] = 1.
        Π[regime][eq[:eq_Ew], ex[:Ew_sh]]   = 1.

        ### E(Rktil)

        # Sticky prices and wages
        Γ0[regime][eq[:eq_ERktil], endo[:Rktil_t]]  = 1.
        Γ1[regime][eq[:eq_ERktil], endo[:ERtil_k_t]] = 1.
        Π[regime][eq[:eq_ERktil], ex[:ERktil_sh]]    = 1.

        ### Wage mark up shock
        if subspec(m) in ["ss52"]
            Γ0[regime][eq[:eq_ϵ_λ_w], endo[:ϵ_λ_w_t]] = 1.
            Γ1[regime][eq[:eq_ϵ_λ_w], endo[:ϵ_λ_w_t]] = m[:η_λ_w]
            if regime == 2
                Γ0[regime][eq[:eq_ϵ_λ_w], endo[:λ_w_t]] = 1. / m[:σ_λ_w_r2]
                Γ1[regime][eq[:eq_ϵ_λ_w], endo[:λ_w_t]] = -m[:ρ_λ_w] / m[:σ_λ_w_r2]
            else
                Γ0[regime][eq[:eq_ϵ_λ_w], endo[:λ_w_t]] = 1. / m[:σ_λ_w]
                Γ1[regime][eq[:eq_ϵ_λ_w], endo[:λ_w_t]] = -m[:ρ_λ_w] / m[:σ_λ_w]
            end
        end

        ### ztil shock
        if subspec(m) in ["ss57"]
            Γ0[regime][eq[:eq_ϵ_ztil], endo[:ϵ_ztil_t]] = 1.
            if regime == 2
                Γ1[regime][eq[:eq_ϵ_ztil], endo[:ztil_t]] = -m[:ρ_ztil] / m[:σ_ztil_r2]
                Γ0[regime][eq[:eq_ϵ_ztil], endo[:z_t]]    = -1 / m[:σ_ztil_r2]
            else
                Γ1[regime][eq[:eq_ϵ_ztil], endo[:ztil_t]] = -m[:ρ_ztil] / m[:σ_ztil]
                Γ0[regime][eq[:eq_ϵ_ztil], endo[:z_t]]    = -1 / m[:σ_ztil]
            end
        end
    end
    return (Γ0[1], Γ1[1], C[1], Ψ[1], Π[1]), (Γ0[2], Γ1[2], C[2], Ψ[2], Π[2])
end
