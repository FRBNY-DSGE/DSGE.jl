# Expresses the equilibrium conditions in canonical form using Γ0, Γ1, C, Ψ, and Π matrices.
# Using the mappings of states/equations to integers defined in m990.jl, coefficients are
# specified in their proper positions.

# Γ0 (n_states x n_states) holds coefficients of current time states.
# Γ1 (n_states x n_states) holds coefficients of lagged states.
# C  (n_states x 1) is a vector of constants
# Ψ  (n_states x n_shocks_exogenous) holds coefficients of iid shocks.
# Π  (n_states x n_states_expectational) holds coefficients of expectational states.
function eqcond(m::SmetsWouters)
    return eqcond(m, 1)
end

function eqcond(m::SmetsWouters, reg::Int)
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
    Γ0[eq[:eq_euler], endo[:c_t]]    = 1.
    Γ0[eq[:eq_euler], endo[:R_t]]    = (1 - m[:h]*exp(-m[:zstar]))/(m[:σ_c]*(1 + m[:h]*exp(-m[:zstar])))
    Γ0[eq[:eq_euler], endo[:b_t]]    = -1.
    Γ0[eq[:eq_euler], endo[:Eπ_t]]   = -(1 - m[:h]*exp(-m[:zstar]))/(m[:σ_c]*(1 + m[:h]*exp(-m[:zstar])))
    Γ0[eq[:eq_euler], endo[:z_t]]    = (m[:h]*exp(-m[:zstar]))/(1 + m[:h]*exp(-m[:zstar]))
    Γ0[eq[:eq_euler], endo[:Ec_t]]   = -1/(1 + m[:h]*exp(-m[:zstar]))
    Γ0[eq[:eq_euler], endo[:ztil_t]] = -( 1/(1-m[:α]) )*(m[:ρ_z]-1)/(1+m[:h]*exp(-m[:zstar]))
    Γ0[eq[:eq_euler], endo[:L_t]]    = -(m[:σ_c] - 1)*m[:wl_c]/(m[:σ_c]*(1 + m[:h]*exp(-m[:zstar])))
    Γ0[eq[:eq_euler], endo[:EL_t]]   = (m[:σ_c] - 1)*m[:wl_c]/(m[:σ_c]*(1 + m[:h]*exp(-m[:zstar])))
    Γ1[eq[:eq_euler], endo[:c_t]]    = (m[:h]*exp(-m[:zstar]))/(1 + m[:h]*exp(-m[:zstar]))

    # Flexible prices and wages
    Γ0[eq[:eq_euler_f], endo[:c_f_t]] = 1.
    Γ0[eq[:eq_euler_f], endo[:r_f_t]] = (1 - m[:h]*exp(-m[:zstar]))/(m[:σ_c]*(1 + m[:h]*exp(-m[:zstar])))
    Γ0[eq[:eq_euler_f], endo[:b_t]]   = -1.
    Γ0[eq[:eq_euler_f], endo[:z_t]]   =   (m[:h]*exp(-m[:zstar]))/(1 + m[:h]*exp(-m[:zstar]))
    Γ0[eq[:eq_euler_f], endo[:Ec_f_t]] = -1/(1 + m[:h]*exp(-m[:zstar]))
    Γ0[eq[:eq_euler_f],endo[:ztil_t]] = -( 1/(1-m[:α]) )*(m[:ρ_z]-1)/(1+m[:h]*exp(-m[:zstar]))
    Γ0[eq[:eq_euler_f], endo[:L_f_t]] = -(m[:σ_c] - 1)*m[:wl_c]/(m[:σ_c]*(1 + m[:h]*exp(-m[:zstar])))
    Γ0[eq[:eq_euler_f], endo[:EL_f_t]] = (m[:σ_c] - 1)*m[:wl_c]/(m[:σ_c]*(1 + m[:h]*exp(-m[:zstar])))
    Γ1[eq[:eq_euler_f], endo[:c_f_t]] = (m[:h]*exp(-m[:zstar]))/(1 + m[:h]*exp(-m[:zstar]))



    ### 2. Investment Euler Equation

    # Sticky prices and wages
    Γ0[eq[:eq_inv], endo[:qk_t]] = -1/(m[:S′′]*exp(2.0*m[:zstar])*(1 + m[:β]*exp((1 - m[:σ_c])*m[:zstar])))
    Γ0[eq[:eq_inv], endo[:i_t]]  = 1.
    Γ0[eq[:eq_inv], endo[:z_t]]  = 1/(1 + m[:β]*exp((1 - m[:σ_c])*m[:zstar]))
    Γ1[eq[:eq_inv], endo[:i_t]]  = 1/(1 + m[:β]*exp((1 - m[:σ_c])*m[:zstar]))
    Γ0[eq[:eq_inv], endo[:Ei_t]]  = -m[:β]*exp((1 - m[:σ_c])*m[:zstar])/(1 + m[:β]*exp((1 - m[:σ_c])*m[:zstar]))
    Γ0[eq[:eq_inv], endo[:ztil_t]] = -( 1/(1-m[:α]) )*(m[:ρ_z]-1)*m[:β]*exp((1-m[:σ_c])*m[:zstar])/(1+m[:β]*exp((1-m[:σ_c])*m[:zstar]))
    Γ0[eq[:eq_inv], endo[:μ_t]] = -1.

    # Flexible prices and wages
    Γ0[eq[:eq_inv_f], endo[:qk_f_t]] = -1/(m[:S′′]*exp(2*m[:zstar])*(1 + m[:β]*exp((1 - m[:σ_c])*m[:zstar])))
    Γ0[eq[:eq_inv_f], endo[:i_f_t]]  = 1.
    Γ0[eq[:eq_inv_f], endo[:z_t]]    = 1/(1 + m[:β]*exp((1 - m[:σ_c])*m[:zstar]))
    Γ1[eq[:eq_inv_f], endo[:i_f_t]]  = 1/(1 + m[:β]*exp((1 - m[:σ_c])*m[:zstar]))
    Γ0[eq[:eq_inv_f], endo[:Ei_f_t]]  = -m[:β]*exp((1 - m[:σ_c])*m[:zstar])/(1 + m[:β]*exp((1 - m[:σ_c])*m[:zstar]))
    Γ0[eq[:eq_inv_f], endo[:ztil_t]] = -( 1/(1-m[:α]) )*(m[:ρ_z]-1)*m[:β]*exp((1-m[:σ_c])*m[:zstar])/(1+m[:β]*exp((1-m[:σ_c])*m[:zstar]))

    Γ0[eq[:eq_inv_f], endo[:μ_t]]   = -1.



    ### 3. Value of Capital

    # Sticky prices and wages

    Γ0[eq[:eq_capval],endo[:Erk_t]] = -m[:rkstar]/(m[:rkstar]+1-m[:δ])
    Γ0[eq[:eq_capval],endo[:Eqk_t]] = -(1-m[:δ])/(m[:rkstar]+1-m[:δ])
    Γ0[eq[:eq_capval],endo[:qk_t]] = 1.
    Γ0[eq[:eq_capval],endo[:R_t]]  = 1.
    Γ0[eq[:eq_capval],endo[:b_t]]  = -(m[:σ_c]*(1+m[:h]*exp(-m[:zstar])))/(1-m[:h]*exp(-m[:zstar]))
    Γ0[eq[:eq_capval],endo[:Eπ_t]] = -1.

    # Flexible prices and wages
    Γ0[eq[:eq_capval_f], endo[:Erk_f_t]] = -m[:rkstar]/(1 + m[:rkstar] - m[:δ])
    Γ0[eq[:eq_capval_f], endo[:Eqk_f_t]] = -(1 - m[:δ])/(1 + m[:rkstar] - m[:δ])
    Γ0[eq[:eq_capval_f], endo[:qk_f_t]] = 1.
    Γ0[eq[:eq_capval_f], endo[:r_f_t]]  = 1.
    Γ0[eq[:eq_capval_f], endo[:b_t]]    = -(m[:σ_c]*(1 + m[:h]*exp(-m[:zstar])))/(1 - m[:h]*exp(-m[:zstar]))



    ### 4. Aggregate Production Function

    # Sticky prices and wages
    Γ0[eq[:eq_output], endo[:y_t]] =  1.
    Γ0[eq[:eq_output], endo[:k_t]] = -m[:Φ]*m[:α]
    Γ0[eq[:eq_output], endo[:L_t]] = -m[:Φ]*(1 - m[:α])
    Γ0[eq[:eq_output], endo[:ztil_t]] = -( m[:Φ]-1 ) / (1-m[:α]);

    # Flexible prices and wages
    Γ0[eq[:eq_output_f], endo[:y_f_t]] =  1.
    Γ0[eq[:eq_output_f], endo[:k_f_t]] = -m[:Φ]*m[:α]
    Γ0[eq[:eq_output_f], endo[:L_f_t]] = -m[:Φ]*(1 - m[:α])
    Γ0[eq[:eq_output_f], endo[:ztil_t]] = -( m[:Φ]-1 ) / (1-m[:α])


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
    Γ0[eq[:eq_capev], endo[:μ_t]]   = -m[:istar]*m[:S′′]*exp(2*m[:zstar])*(1 + m[:β]*exp((1 - m[:σ_c])*m[:zstar]))/m[:kbarstar]

    # Flexible prices and wages
    Γ0[eq[:eq_capev_f], endo[:kbar_f_t]] = 1.
    Γ1[eq[:eq_capev_f], endo[:kbar_f_t]] = 1 - m[:istar]/m[:kbarstar]
    Γ0[eq[:eq_capev_f], endo[:z_t]]      = 1 - m[:istar]/m[:kbarstar]
    Γ0[eq[:eq_capev_f], endo[:i_f_t]]    = -m[:istar]/m[:kbarstar]
    Γ0[eq[:eq_capev_f], endo[:μ_t]]     = -m[:istar]*m[:S′′]*exp(2*m[:zstar])*(1 + m[:β]*exp((1 - m[:σ_c])*m[:zstar]))/m[:kbarstar]



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
       Γ0[eq[:eq_phlps], endo[:π_t]] = 1.
    Γ0[eq[:eq_phlps], endo[:mc_t]] =  -((1 - m[:ζ_p]*m[:β]*exp((1 - m[:σ_c])*m[:zstar]))*(1 - m[:ζ_p]))/(m[:ζ_p]*((m[:Φ]- 1)*m[:ϵ_p] + 1))/(1 + m[:ι_p]*m[:β]*exp((1 - m[:σ_c])*m[:zstar]))
    Γ1[eq[:eq_phlps], endo[:π_t]] = m[:ι_p]/(1 + m[:ι_p]*m[:β]*exp((1 - m[:σ_c])*m[:zstar]))
    Γ0[eq[:eq_phlps], endo[:Eπ_t]] = -m[:β]*exp((1 - m[:σ_c])*m[:zstar])/(1 + m[:ι_p]*m[:β]*exp((1 - m[:σ_c])*m[:zstar]))

    Γ0[eq[:eq_phlps], endo[:λ_f_t]] = -(1+m[:ι_p]*m[:β]*exp((1-m[:σ_c])*m[:zstar]))/(1+m[:ι_p]*m[:β]*exp((1-m[:σ_c])*m[:zstar]))

    # Flexible prices and wages not necessary

    ### 10. Rental Rate of Capital

    # Sticky prices and wages
    Γ0[eq[:eq_caprnt], endo[:rk_t]] = 1.
    Γ0[eq[:eq_caprnt], endo[:k_t]]  = 1.
    Γ0[eq[:eq_caprnt], endo[:L_t]]  = -1.
    Γ0[eq[:eq_caprnt], endo[:w_t]]  = -1.

    # Flexible prices and wages
    Γ0[eq[:eq_caprnt_f], endo[:rk_f_t]] = 1.
    Γ0[eq[:eq_caprnt_f], endo[:k_f_t]] = 1.
    Γ0[eq[:eq_caprnt_f], endo[:L_f_t]] = -1.
    Γ0[eq[:eq_caprnt_f], endo[:w_f_t]] = -1.



    ### 11. Marginal Substitution

    # Sticky prices and wages
    Γ0[eq[:eq_msub], endo[:μ_ω_t]] = 1.
    Γ0[eq[:eq_msub], endo[:L_t]]   = m[:ν_l]
    Γ0[eq[:eq_msub], endo[:c_t]]   = 1/(1 - m[:h]*exp(-m[:zstar]))
    Γ1[eq[:eq_msub], endo[:c_t]]   = m[:h]*exp(-m[:zstar])/(1 - m[:h]*exp(-m[:zstar]))
    Γ0[eq[:eq_msub], endo[:z_t]]   = m[:h]*exp(-m[:zstar]) /(1 - m[:h]*exp(-m[:zstar]))
    Γ0[eq[:eq_msub], endo[:w_t]]   = -1.

    # Flexible prices and wages
    Γ0[eq[:eq_msub_f], endo[:w_f_t]] = -1.
    Γ0[eq[:eq_msub_f], endo[:L_f_t]] = m[:ν_l]
    Γ0[eq[:eq_msub_f], endo[:c_f_t]] = 1/(1 - m[:h]*exp(-m[:zstar]))
    Γ1[eq[:eq_msub_f], endo[:c_f_t]] = m[:h]*exp(-m[:zstar])/(1 - m[:h]*exp(-m[:zstar]))
    Γ0[eq[:eq_msub_f], endo[:z_t]]   = m[:h]*exp(-m[:zstar])/(1 - m[:h]*exp(-m[:zstar]))


    ### 12. Evolution of Wages

    # Sticky prices and wages
    Γ0[eq[:eq_wage], endo[:w_t]]   = 1
    Γ0[eq[:eq_wage], endo[:μ_ω_t]] = (1 - m[:ζ_w]*m[:β]*exp((1 - m[:σ_c])*m[:zstar]))*
        (1 - m[:ζ_w])/(m[:ζ_w]*((m[:λ_w] - 1)*m[:ϵ_w] + 1))/(1 + m[:β]*exp((1 - m[:σ_c])*m[:zstar]))
    Γ0[eq[:eq_wage], endo[:π_t]]  = (1 + m[:ι_w]*m[:β]*exp((1 - m[:σ_c])*m[:zstar]))/(1 + m[:β]*exp((1 - m[:σ_c])*m[:zstar]))
    Γ1[eq[:eq_wage], endo[:w_t]]   = 1/(1 + m[:β]*exp((1 - m[:σ_c])*m[:zstar]))
    Γ0[eq[:eq_wage], endo[:z_t]]   = 1/(1 + m[:β]*exp((1 - m[:σ_c])*m[:zstar]))
    Γ1[eq[:eq_wage], endo[:π_t]]  = m[:ι_w]/(1 + m[:β]*exp((1 - m[:σ_c])*m[:zstar]))
    Γ0[eq[:eq_wage], endo[:Ew_t]]   = -m[:β]*exp((1 - m[:σ_c])*m[:zstar])/(1 + m[:β]*exp((1 - m[:σ_c])*m[:zstar]))
    Γ0[eq[:eq_wage], endo[:ztil_t]] = -( 1/(1-m[:α]) )*(m[:ρ_z]-1)*m[:β]*exp((1-m[:σ_c])*m[:zstar])*1/(1+m[:β]*exp((1-m[:σ_c])*m[:zstar]))
    Γ0[eq[:eq_wage], endo[:Eπ_t]]  = -m[:β]*exp((1 - m[:σ_c])*m[:zstar])/(1 + m[:β]*exp((1 - m[:σ_c])*m[:zstar]))
    Γ0[eq[:eq_wage], endo[:λ_w_t]] = -1.

    # Flexible prices and wages not necessary



    ### 13. Monetary Policy Rule

    # Sticky prices and wages
    Γ0[eq[:eq_mp], endo[:R_t]]    = 1.
    Γ1[eq[:eq_mp], endo[:R_t]]    = m[:ρ]
    Γ0[eq[:eq_mp], endo[:π_t]]   = -(1 - m[:ρ])*m[:ψ1]
    Γ0[eq[:eq_mp], endo[:y_t]]    = -(1 - m[:ρ])*m[:ψ2] - m[:ψ3]
    Γ0[eq[:eq_mp], endo[:y_f_t]]  = (1 - m[:ρ])*m[:ψ2] + m[:ψ3]
    Γ1[eq[:eq_mp], endo[:y_t]]    = -m[:ψ3]
    Γ1[eq[:eq_mp], endo[:y_f_t]]  = m[:ψ3]
    Γ0[eq[:eq_mp], endo[:rm_t]]   = -1.

    # Flexible prices and wages not necessary



    ### 14. Resource Constraint

    # Sticky prices and wages
    Γ0[eq[:eq_res], endo[:y_t]] = 1.
    Γ0[eq[:eq_res], endo[:g_t]] = -m[:g_star]
    if m[:ρ_z] < 1
        Γ0[eq[:eq_res], endo[:ztil_t]] = m[:g_star] * (1/(1-m[:α]))
    end
    Γ0[eq[:eq_res], endo[:c_t]] = -m[:cstar]/m[:ystar]
    Γ0[eq[:eq_res], endo[:i_t]] = -m[:istar]/m[:ystar]
    Γ0[eq[:eq_res], endo[:u_t]] = -m[:rkstar]*m[:kstar]/m[:ystar]

    # Flexible prices and wages
    Γ0[eq[:eq_res_f], endo[:y_f_t]] = 1.
    Γ0[eq[:eq_res_f], endo[:g_t]]   = -m[:g_star]
    if m[:ρ_z]<1
        Γ0[eq[:eq_res_f], endo[:ztil_t]] = m[:g_star] * (1/(1-m[:α]))
    end
    Γ0[eq[:eq_res_f], endo[:c_f_t]] = -m[:cstar]/m[:ystar]
    Γ0[eq[:eq_res_f], endo[:i_f_t]] = -m[:istar]/m[:ystar]
    Γ0[eq[:eq_res_f], endo[:u_f_t]] = -m[:rkstar]*m[:kstar]/m[:ystar]


    ### EXOGENOUS SHOCKS ###

    # Neutral technology
    Γ0[eq[:eq_z], endo[:z_t]]    = 1.
    Γ1[eq[:eq_z], endo[:ztil_t]] = (m[:ρ_z] - 1)/(1 - m[:α])
    Ψ[eq[:eq_z], exo[:z_sh]]     = 1/(1 - m[:α])

    Γ0[eq[:eq_ztil], endo[:ztil_t]] = 1.
    Γ1[eq[:eq_ztil], endo[:ztil_t]] = m[:ρ_z]
    Ψ[eq[:eq_ztil], exo[:z_sh]]     = 1.

    # Government spending
    Γ0[eq[:eq_g], endo[:g_t]] = 1.
    Γ1[eq[:eq_g], endo[:g_t]] = m[:ρ_g]
    Ψ[eq[:eq_g], exo[:g_sh]]  = 1.
    Ψ[eq[:eq_g], exo[:z_sh]]  = m[:η_gz]

    # Asset shock
    Γ0[eq[:eq_b], endo[:b_t]] = 1.
    Γ1[eq[:eq_b], endo[:b_t]] = m[:ρ_b]
    Ψ[eq[:eq_b], exo[:b_sh]]  = 1.

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

    # Wage mark-up shock
    Γ0[eq[:eq_λ_w], endo[:λ_w_t]]  = 1.
    Γ1[eq[:eq_λ_w], endo[:λ_w_t]]  = m[:ρ_λ_w]
    Γ1[eq[:eq_λ_w], endo[:λ_w_t1]] = -m[:η_λ_w]
    Ψ[eq[:eq_λ_w], exo[:λ_w_sh]]   = 1.

    Γ0[eq[:eq_λ_w1], endo[:λ_w_t1]] = 1.
    Ψ[eq[:eq_λ_w1], exo[:λ_w_sh]]   = 1.

    # Monetary policy shock
    Γ0[eq[:eq_rm], endo[:rm_t]] = 1.
    Γ1[eq[:eq_rm], endo[:rm_t]] = m[:ρ_rm]
    Ψ[eq[:eq_rm], exo[:rm_sh]]  = 1.


    # Anticipated policy shocks
    if n_anticipated_shocks(m) > 0

        # This section adds the anticipated shocks. There is one state for all the
        # anticipated shocks that will hit in a given period (i.e. rm_tl2 holds those that
        # will hit in two periods), and the equations are set up so that rm_tl2 last period
        # will feed into rm_tl1 this period (and so on for other numbers), and last period's
        # rm_tl1 will feed into the rm_t process (and affect the Taylor Rule this period).

        Γ1[eq[:eq_rm], endo[:rm_tl1]]   = 1.
        Γ0[eq[:eq_rml1], endo[:rm_tl1]] = 1.
        Ψ[eq[:eq_rml1], exo[:rm_shl1]]  = 1.

        if n_anticipated_shocks(m) > 1
            for i = 2:n_anticipated_shocks(m)
                Γ1[eq[Symbol("eq_rml$(i-1)")], endo[Symbol("rm_tl$i")]] = 1.
                Γ0[eq[Symbol("eq_rml$i")], endo[Symbol("rm_tl$i")]] = 1.
                Ψ[eq[Symbol("eq_rml$i")], exo[Symbol("rm_shl$i")]] = 1.
            end
        end
    end



    ### EXPECTATION ERRORS ###

    ### E(c)

    # Sticky prices and wages
    Γ0[eq[:eq_Ec], endo[:c_t]]         = 1.
    Γ1[eq[:eq_Ec], endo[:Ec_t]]        = 1.
    Π[eq[:eq_Ec], ex[:Ec_sh]]          = 1.

    # Flexible prices and wages
    Γ0[eq[:eq_Ec_f], endo[:c_f_t]]     = 1.
    Γ1[eq[:eq_Ec_f], endo[:Ec_f_t]]    = 1.
    Π[eq[:eq_Ec_f], ex[:Ec_f_sh]]      = 1.



    ### E(q)

    # Sticky prices and wages
    Γ0[eq[:eq_Eqk], endo[:qk_t]]       = 1.
    Γ1[eq[:eq_Eqk], endo[:Eqk_t]]      = 1.
    Π[eq[:eq_Eqk], ex[:Eqk_sh]]        = 1.

    # Flexible prices and wages
    Γ0[eq[:eq_Eqk_f], endo[:qk_f_t]]   = 1.
    Γ1[eq[:eq_Eqk_f], endo[:Eqk_f_t]]  = 1.
    Π[eq[:eq_Eqk_f], ex[:Eqk_f_sh]]    = 1.

    ### E(i)

    # Sticky prices and wages
    Γ0[eq[:eq_Ei], endo[:i_t]]         = 1.
    Γ1[eq[:eq_Ei], endo[:Ei_t]]         = 1.
    Π[eq[:eq_Ei], ex[:Ei_sh]]          = 1.

    # Flexible prices and wages
    Γ0[eq[:eq_Ei_f], endo[:i_f_t]]     = 1.
    Γ1[eq[:eq_Ei_f], endo[:Ei_f_t]]     = 1.
    Π[eq[:eq_Ei_f], ex[:Ei_f_sh]]      = 1.



    ### E(π)

    # Sticky prices and wages
    Γ0[eq[:eq_Eπ], endo[:π_t]]       = 1.
    Γ1[eq[:eq_Eπ], endo[:Eπ_t]]       = 1.
    Π[eq[:eq_Eπ], ex[:Eπ_sh]]        = 1.



    ### E(l)

    # Sticky prices and wages
    Γ0[eq[:eq_EL], endo[:L_t]]         = 1.
    Γ1[eq[:eq_EL], endo[:EL_t]]         = 1.
    Π[eq[:eq_EL], ex[:EL_sh]]          = 1.

    # Flexible prices and wages
    Γ0[eq[:eq_EL_f], endo[:L_f_t]]     = 1.
    Γ1[eq[:eq_EL_f], endo[:EL_f_t]]     = 1.
    Π[eq[:eq_EL_f], ex[:EL_f_sh]]      = 1.



    ### E(rk)

    # Sticky prices and wages
    Γ0[eq[:eq_Erk], endo[:rk_t]]       = 1.
    Γ1[eq[:eq_Erk], endo[:Erk_t]]       = 1.
    Π[eq[:eq_Erk], ex[:Erk_sh]]        = 1.

    # Flexible prices and wages
    Γ0[eq[:eq_Erk_f], endo[:rk_f_t]]   = 1.
    Γ1[eq[:eq_Erk_f], endo[:Erk_f_t]]   = 1.
    Π[eq[:eq_Erk_f], ex[:Erk_f_sh]]    = 1.



    ### E(w)

    # Sticky prices and wages
    Γ0[eq[:eq_Ew], endo[:w_t]]         = 1.
    Γ1[eq[:eq_Ew], endo[:Ew_t]]         = 1.
    Π[eq[:eq_Ew], ex[:Ew_sh]]          = 1.

    for para in m.parameters
        if !isempty(para.regimes)
            ModelConstructors.toggle_regime!(para, 1)
        end
    end

    return Γ0, Γ1, C, Ψ, Π
end
