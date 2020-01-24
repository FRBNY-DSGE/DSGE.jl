#g Expresses the equilibrium conditions in canonical form using Γ0, Γ1, C, Ψ, and Π matrices.
# Using the mappings of states/equations to integers defined in m990.jl, coefficients are
# specified in their proper positions.

# Γ0 (n_states x n_states) holds coefficients of current time states.
# Γ1 (n_states x n_states) holds coefficients of lagged states.
# C  (n_states x 1) is a vector of constants
# Ψ  (n_states x n_shocks_exogenous) holds coefficients of iid shocks.
# Π  (n_states x n_states_expectational) holds coefficients of expectational states.

function eqcond_regimes(m::SmetsWoutersOrig)
    Γ0 = Vector{Array}(undef, 2)
    Γ1 = Vector{Array}(undef, 2)
    C = Vector{Array}(undef, 2)
    Ψ = Vector{Array}(undef, 2)
    Π = Vector{Array}(undef, 2)

  #  if subspec(m) != "ss51"
        for regime = 1:2
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

            #TRYING TO MATCH TO FORTRAN
            #m[:γ] = m[:γ]/100 + 1
            #m[:β] = 1/(m[:β]/100+1)
            #m[:π_star] = m[:π_star]/100+1

            # Sticky prices and wages
            Γ0[regime][eq[:eq_euler], endo[:c_t]]    = 1.
            Γ0[regime][eq[:eq_euler], endo[:Ec_t]]   = -1/(1 + m[:h]/m[:γ])
            Γ0[regime][eq[:eq_euler], endo[:L_t]]    = -(m[:σ_c] - 1)*m[:wl_c]/(m[:σ_c]*(1 + m[:h]/m[:γ]))
            Γ0[regime][eq[:eq_euler], endo[:EL_t]]   = (m[:σ_c] - 1)*m[:wl_c]/(m[:σ_c]*(1 + m[:h]/m[:γ]))
            Γ0[regime][eq[:eq_euler], endo[:R_t]]    = (1 - m[:h]/m[:γ])/(m[:σ_c]*(1 + m[:h]/m[:γ]))
            Γ0[regime][eq[:eq_euler], endo[:Eπ_t]]   = -(1 - m[:h]/m[:γ])/(m[:σ_c]*(1 + m[:h]/m[:γ]))
            Γ0[regime][eq[:eq_euler], endo[:b_t]]    = -1.
            Γ1[regime][eq[:eq_euler], endo[:c_t]]    = (m[:h]/m[:γ])/(1 + m[:h]/m[:γ])

            # Flexible prices and wages
            Γ0[regime][eq[:eq_euler_f], endo[:c_f_t]]  = 1.
            Γ0[regime][eq[:eq_euler_f], endo[:Ec_f_t]] = -1/(1 + m[:h]/m[:γ])
            Γ0[regime][eq[:eq_euler_f], endo[:L_f_t]]  = -(m[:σ_c] - 1)*m[:wl_c]/(m[:σ_c]*(1 + m[:h]/m[:γ]))
            Γ0[regime][eq[:eq_euler_f], endo[:EL_f_t]] = (m[:σ_c] - 1)*m[:wl_c]/(m[:σ_c]*(1 + m[:h]/m[:γ]))
            Γ0[regime][eq[:eq_euler_f], endo[:r_f_t]]  = (1 - m[:h]/m[:γ])/(m[:σ_c]*(1 + m[:h]/m[:γ]))
            Γ0[regime][eq[:eq_euler_f], endo[:b_t]]    = -1.
            Γ1[regime][eq[:eq_euler_f], endo[:c_f_t]]  = (m[:h]/m[:γ])/(1 + m[:h]/m[:γ])


            ### 2. Investment Euler Equation

            # Sticky prices and wages
            Γ0[regime][eq[:eq_inv], endo[:i_t]]  = 1.
            Γ0[regime][eq[:eq_inv], endo[:Ei_t]] = -m[:β]*m[:γ]^(1 - m[:σ_c])/(1 + m[:β]*m[:γ]^(1 - m[:σ_c]))
            Γ0[regime][eq[:eq_inv], endo[:qk_t]] = -1 / (m[:S′′] * m[:γ]^2 * (1 + m[:β]*m[:γ]^(1 - m[:σ_c])))
            Γ0[regime][eq[:eq_inv], endo[:μ_t]]  = -1.
            Γ1[regime][eq[:eq_inv], endo[:i_t]]  = 1/(1 + m[:β]*m[:γ]^(1 - m[:σ_c]))

            # Flexible prices and wages
            Γ0[regime][eq[:eq_inv_f], endo[:i_f_t]]  = 1.
            Γ0[regime][eq[:eq_inv_f], endo[:Ei_f_t]] = -m[:β]*m[:γ]^(1 - m[:σ_c])/(1 + m[:β]*m[:γ]^(1 - m[:σ_c]))
            Γ0[regime][eq[:eq_inv_f], endo[:qk_f_t]] = -1 / (m[:S′′] * m[:γ]^2 * (1 + m[:β]*m[:γ]^(1 - m[:σ_c])))
            Γ0[regime][eq[:eq_inv_f], endo[:μ_t]]    = -1.
            Γ1[regime][eq[:eq_inv_f], endo[:i_f_t]]  = 1/(1 + m[:β]*m[:γ]^(1 - m[:σ_c]))


            ### 3. Value of Capital

            # Sticky prices and wages

            Γ0[regime][eq[:eq_capval],endo[:qk_t]]  = 1.
            Γ0[regime][eq[:eq_capval],endo[:Eqk_t]] = -m[:β] * m[:γ]^(-m[:σ_c]) * (1 - m[:δ])
            Γ0[regime][eq[:eq_capval],endo[:Erk_t]] = -(1 - m[:β] * m[:γ]^(-m[:σ_c]) * (1 - m[:δ]))
            Γ0[regime][eq[:eq_capval],endo[:R_t]]   = 1.
            Γ0[regime][eq[:eq_capval],endo[:Eπ_t]]  = -1.
            Γ0[regime][eq[:eq_capval],endo[:b_t]]   = -1 / ((1-m[:h]/m[:γ]) / ((1+m[:h]/m[:γ])* m[:σ_c]))

            # Flexible prices and wages
            Γ0[regime][eq[:eq_capval_f], endo[:qk_f_t]]  = 1.
            Γ0[regime][eq[:eq_capval_f], endo[:Eqk_f_t]] = -m[:β] * m[:γ]^(-m[:σ_c]) * (1 - m[:δ])
            Γ0[regime][eq[:eq_capval_f], endo[:Erk_f_t]] = -(1 - m[:β] * m[:γ]^(-m[:σ_c]) * (1 - m[:δ]))
            Γ0[regime][eq[:eq_capval_f], endo[:r_f_t]]   = 1.
            Γ0[regime][eq[:eq_capval_f], endo[:b_t]]     = -1 / ((1-m[:h]/m[:γ]) / ((1+m[:h]/m[:γ])* m[:σ_c]))


            ### 4. Aggregate Production Function

            # Sticky prices and wages
            Γ0[regime][eq[:eq_output], endo[:y_t]] =  1.
            Γ0[regime][eq[:eq_output], endo[:k_t]] = -m[:Φ]*m[:α]
            Γ0[regime][eq[:eq_output], endo[:L_t]] = -m[:Φ]*(1 - m[:α])
            Γ0[regime][eq[:eq_output], endo[:z_t]] = -m[:Φ]

            # Flexible prices and wages
            Γ0[regime][eq[:eq_output_f], endo[:y_f_t]] =  1.
            Γ0[regime][eq[:eq_output_f], endo[:k_f_t]] = -m[:Φ]*m[:α]
            Γ0[regime][eq[:eq_output_f], endo[:L_f_t]] = -m[:Φ]*(1 - m[:α])
            Γ0[regime][eq[:eq_output_f], endo[:z_t]]   = -m[:Φ]


    ### 5. Capital Utilization

    # Sticky prices and wages
        Γ0[regime][eq[:eq_caputl], endo[:k_t]]    =  1.
        Γ1[regime][eq[:eq_caputl], endo[:kbar_t]] =  1.
        Γ0[regime][eq[:eq_caputl], endo[:u_t]]    = -1.

    # Flexible prices and wages
        Γ0[regime][eq[:eq_caputl_f], endo[:k_f_t]]    =  1.
        Γ1[regime][eq[:eq_caputl_f], endo[:kbar_f_t]] =  1.
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
        Γ0[regime][eq[:eq_capev], endo[:μ_t]]    = -(1 - (1 - m[:δ])/m[:γ]) * m[:γ]^2 * m[:S′′] * (1 + m[:β]*m[:γ]^(1 - m[:σ_c]))
        Γ0[regime][eq[:eq_capev], endo[:i_t]]    = -(1 - (1 - m[:δ])/m[:γ])
        Γ1[regime][eq[:eq_capev], endo[:kbar_t]] = (1 - m[:δ])/m[:γ]

    # Flexible prices and wages
        Γ0[regime][eq[:eq_capev_f], endo[:kbar_f_t]] = 1.
        Γ0[regime][eq[:eq_capev_f], endo[:μ_t]]      = -(1 - (1 - m[:δ])/m[:γ]) * m[:γ]^2 * m[:S′′] * (1 + m[:β]*m[:γ]^(1 - m[:σ_c]))
        Γ0[regime][eq[:eq_capev_f], endo[:i_f_t]]    = -(1 - (1 - m[:δ])/m[:γ])
        Γ1[regime][eq[:eq_capev_f], endo[:kbar_f_t]] = (1 - m[:δ])/m[:γ]


    ### 8. Price Markup

    # Sticky prices and wages
        Γ0[regime][eq[:eq_mkupp], endo[:mc_t]] =  1.
        Γ0[regime][eq[:eq_mkupp], endo[:k_t]]  = -m[:α]
        Γ0[regime][eq[:eq_mkupp], endo[:L_t]]  =  m[:α]
        Γ0[regime][eq[:eq_mkupp], endo[:z_t]]  = -1.
        Γ0[regime][eq[:eq_mkupp], endo[:w_t]]  =  1.

    # Flexible prices and wages
        Γ0[regime][eq[:eq_mkupp_f], endo[:k_f_t]] = -m[:α]
        Γ0[regime][eq[:eq_mkupp_f], endo[:L_f_t]] =  m[:α]
        Γ0[regime][eq[:eq_mkupp_f], endo[:z_t]]   = -1.
        Γ0[regime][eq[:eq_mkupp_f], endo[:w_f_t]] =  1.


    ### 9. Phillips Curve

    # Sticky prices and wages
        Γ0[regime][eq[:eq_phlps], endo[:π_t]]   = 1.
        Γ0[regime][eq[:eq_phlps], endo[:Eπ_t]]  = -m[:β] * m[:γ]^(1 - m[:σ_c]) / (1 + m[:β]*m[:γ]^(1 - m[:σ_c])*m[:ι_p])
        if subspec(m) in ["ss21", "ss22", "ss25", "ss26", "ss28", "ss29", "ss41", "ss42"] && regime == 2
            Γ0[regime][eq[:eq_phlps], endo[:mc_t]]  = ((1 - m[:ζ_p_r2])*(1 - m[:ζ_p_r2]*m[:β]*m[:γ]^(1 - m[:σ_c]))) /
                (m[:ζ_p_r2]*((m[:Φ]- 1)*m[:ϵ_p] + 1)*(1 + m[:β]*m[:γ]^(1 - m[:σ_c])*m[:ι_p]))
        else
            Γ0[regime][eq[:eq_phlps], endo[:mc_t]]  = ((1 - m[:ζ_p])*(1 - m[:ζ_p]*m[:β]*m[:γ]^(1 - m[:σ_c]))) /
                (m[:ζ_p]*((m[:Φ]- 1)*m[:ϵ_p] + 1)*(1 + m[:β]*m[:γ]^(1 - m[:σ_c])*m[:ι_p]))
        end

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
            κnum = ((1 - m[:ζ_p]*m[:β]*exp((1 - m[:σ_c])))* #m[:zstar]))*
                    (1 - m[:ζ_p]))/(m[:ζ_p]*((m[:Φ]- 1)*m[:ϵ_p] + 1))/(1 + m[:ι_p]*m[:β]*exp((1 - m[:σ_c]))) #*m[:zstar]))         # kappa numerator
            fix_ζ_p = m[:ζ_p2]
            κden = ((1 - fix_ζ_p*m[:β]*exp((1 - m[:σ_c])))* #m[:zstar]))*
                    (1 - fix_ζ_p))/(fix_ζ_p*((m[:Φ]- 1)*m[:ϵ_p] + 1))/(1 + m[:ι_p]*m[:β]*exp((1 - m[:σ_c]))) #*m[:zstar]))         # kappa denominator
            Γ0[regime][eq[:eq_phlps], endo[:λ_f_t]] = -κnum / κden
       #= elseif subspec(m) in ["ss21", "ss22", "ss25", "ss26", "ss28", "ss29", "ss41", "ss42"]
            if regime == 1
                Γ0[regime][eq[:eq_phlps], endo[:λ_f_t]] = -1.
            elseif regime == 2
                κnum = ((1 - m[:ζ_p_r2])*(1 - m[:ζ_p_r2]*m[:β]*m[:γ]^(1 - m[:σ_c]))) /
                    (m[:ζ_p_r2]*((m[:Φ]- 1)*m[:ϵ_p] + 1)*(1 + m[:β]*m[:γ]^(1 - m[:σ_c])*m[:ι_p]))
                fix_ζ_p = m[:ζ_p]
                κden = ((1 - fix_ζ_p)*(1 - fix_ζ_p*m[:β]*m[:γ]^(1 - m[:σ_c]))) /
                    (fix_ζ_p*((m[:Φ]- 1)*m[:ϵ_p] + 1)*(1 + m[:β]*m[:γ]^(1 - m[:σ_c])*m[:ι_p]))
                Γ0[regime][eq[:eq_phlps], endo[:λ_f_t]] = -κnum / κden
            end =#
        else
            Γ0[regime][eq[:eq_phlps], endo[:λ_f_t]] = -1.
        end
        Γ1[regime][eq[:eq_phlps], endo[:π_t]]   = m[:ι_p] / (1 + m[:β]*m[:γ]^(1 - m[:σ_c])*m[:ι_p])

    # Flexible prices and wages not necessary

    ### 10. Rental Rate of Capital

    # Sticky prices and wages
            Γ0[regime][eq[:eq_caprnt], endo[:rk_t]] =  1.
            Γ0[regime][eq[:eq_caprnt], endo[:k_t]]  =  1.
            Γ0[regime][eq[:eq_caprnt], endo[:L_t]]  = -1.
            Γ0[regime][eq[:eq_caprnt], endo[:w_t]]  = -1.

    # Flexible prices and wages
            Γ0[regime][eq[:eq_caprnt_f], endo[:rk_f_t]] =  1.
            Γ0[regime][eq[:eq_caprnt_f], endo[:k_f_t]]  =  1.
            Γ0[regime][eq[:eq_caprnt_f], endo[:L_f_t]]  = -1.
            Γ0[regime][eq[:eq_caprnt_f], endo[:w_f_t]]  = -1.


    ### 11. Wage Markup (Marginal Substitution)

    # Sticky prices and wages
            Γ0[regime][eq[:eq_msub], endo[:μ_ω_t]] =  1.
            Γ0[regime][eq[:eq_msub], endo[:w_t]]   = -1.
            Γ0[regime][eq[:eq_msub], endo[:L_t]]   = m[:ν_l]
            Γ0[regime][eq[:eq_msub], endo[:c_t]]   = 1/(1 - m[:h]/m[:γ])
            Γ1[regime][eq[:eq_msub], endo[:c_t]]   = (m[:h]/m[:γ]) / (1 - m[:h]/m[:γ])

    # Flexible prices and wages
            Γ0[regime][eq[:eq_msub_f], endo[:w_f_t]] = -1.
            Γ0[regime][eq[:eq_msub_f], endo[:L_f_t]] = m[:ν_l]
            Γ0[regime][eq[:eq_msub_f], endo[:c_f_t]] = 1/(1 - m[:h]/m[:γ])
            Γ1[regime][eq[:eq_msub_f], endo[:c_f_t]] = (m[:h]/m[:γ]) / (1 - m[:h]/m[:γ])


    ### 12. Evolution of Wages

    # Sticky prices and wages
            Γ0[regime][eq[:eq_wage], endo[:w_t]]    = 1
            Γ0[regime][eq[:eq_wage], endo[:Ew_t]]   = -(1 - 1/(1 + m[:β]*m[:γ]^(1 - m[:σ_c])))
            Γ0[regime][eq[:eq_wage], endo[:Eπ_t]]   = -(1 - 1/(1 + m[:β]*m[:γ]^(1 - m[:σ_c])))
            Γ0[regime][eq[:eq_wage], endo[:π_t]]    = (1 + m[:β]*m[:γ]^(1 - m[:σ_c])*m[:ι_w]) / (1 + m[:β]*m[:γ]^(1 - m[:σ_c]))
            Γ0[regime][eq[:eq_wage], endo[:μ_ω_t]]  = ((1 - m[:β]*m[:γ]^(1 - m[:σ_c])*m[:ζ_w]) * (1 - m[:ζ_w])) /
        ((1 + m[:β]*m[:γ]^(1 - m[:σ_c])) * (m[:ζ_w]*((m[:λ_w] - 1)*m[:ϵ_w] + 1)))
            Γ0[regime][eq[:eq_wage], endo[:λ_w_t]]  = -1.
            Γ1[regime][eq[:eq_wage], endo[:w_t]]    = 1/(1 + m[:β]*m[:γ]^(1 - m[:σ_c]))
            Γ1[regime][eq[:eq_wage], endo[:π_t]]    = m[:ι_w]/(1 + m[:β]*m[:γ]^(1 - m[:σ_c]))

    # Flexible prices and wages not necessary


    ### 13. Monetary Policy Rule

    # Sticky prices and wages
            Γ0[regime][eq[:eq_mp], endo[:R_t]]   = 1.
        if subspec(m) in ["ss23", "ss24", "ss25", "ss26", "ss28", "ss29", "ss43", "ss44"] && regime == 2
                Γ1[regime][eq[:eq_mp], endo[:R_t]]   = m[:ρ_r2]
                Γ0[regime][eq[:eq_mp], endo[:π_t]]   = -(1 - m[:ρ_r2])*m[:ψ1_r2]
                Γ0[regime][eq[:eq_mp], endo[:y_t]]   = -(1 - m[:ρ_r2])*m[:ψ2_r2] - m[:ψ3_r2]
                Γ1[regime][eq[:eq_mp], endo[:y_t]]   = -m[:ψ3_r2]
                Γ0[regime][eq[:eq_mp], endo[:y_f_t]] = (1 - m[:ρ_r2])*m[:ψ2_r2] + m[:ψ3_r2]
                Γ1[regime][eq[:eq_mp], endo[:y_f_t]] = m[:ψ3_r2]
                Γ0[regime][eq[:eq_mp], endo[:rm_t]]  = -1.
        else
                Γ1[regime][eq[:eq_mp], endo[:R_t]]   = m[:ρ]
                Γ0[regime][eq[:eq_mp], endo[:π_t]]   = -(1 - m[:ρ])*m[:ψ1]
                Γ0[regime][eq[:eq_mp], endo[:y_t]]   = -(1 - m[:ρ])*m[:ψ2] - m[:ψ3]
                Γ1[regime][eq[:eq_mp], endo[:y_t]]   = -m[:ψ3]
                Γ0[regime][eq[:eq_mp], endo[:y_f_t]] = (1 - m[:ρ])*m[:ψ2] + m[:ψ3]
                Γ1[regime][eq[:eq_mp], endo[:y_f_t]] = m[:ψ3]
                Γ0[regime][eq[:eq_mp], endo[:rm_t]]  = -1.
        end
    # Flexible prices and wages not necessary


    ### 14. Resource Constraint

    # Sticky prices and wages
            Γ0[regime][eq[:eq_res], endo[:y_t]] =  1.
            Γ0[regime][eq[:eq_res], endo[:c_t]] = -m[:c_y]
            Γ0[regime][eq[:eq_res], endo[:i_t]] = -m[:i_y]
            Γ0[regime][eq[:eq_res], endo[:u_t]] = -m[:u_y]
            Γ0[regime][eq[:eq_res], endo[:g_t]] = -1.

    # Flexible prices and wages
            Γ0[regime][eq[:eq_res_f], endo[:y_f_t]] = 1.
            Γ0[regime][eq[:eq_res_f], endo[:c_f_t]] = -m[:c_y]
            Γ0[regime][eq[:eq_res_f], endo[:i_f_t]] = -m[:i_y]
            Γ0[regime][eq[:eq_res_f], endo[:u_f_t]] = -m[:u_y]
            Γ0[regime][eq[:eq_res_f], endo[:g_t]]   = -1.


    ### EXOGENOUS PROCESSES ###

    # Government spending
            Γ0[regime][eq[:eq_g], endo[:g_t]] = 1.
            Γ1[regime][eq[:eq_g], endo[:g_t]] = m[:ρ_g]
        Ψ[regime][eq[:eq_g], exo[:g_sh]]  = 1.
        Ψ[regime][eq[:eq_g], exo[:z_sh]]  = m[:η_gz]

    # Asset shock
            Γ0[regime][eq[:eq_b], endo[:b_t]] = 1.
            Γ1[regime][eq[:eq_b], endo[:b_t]] = m[:ρ_b]
        Ψ[regime][eq[:eq_b], exo[:b_sh]]  = 1.

    # Investment-specific technology
            Γ0[regime][eq[:eq_μ], endo[:μ_t]] = 1.
            Γ1[regime][eq[:eq_μ], endo[:μ_t]] = m[:ρ_μ]
        Ψ[regime][eq[:eq_μ], exo[:μ_sh]]  = 1.

    # Neutral technology
            Γ0[regime][eq[:eq_z], endo[:z_t]] = 1.
            Γ1[regime][eq[:eq_z], endo[:z_t]] = m[:ρ_z]
        Ψ[regime][eq[:eq_z], exo[:z_sh]]  = 1.

    # Price mark-up shock
            Γ0[regime][eq[:eq_λ_f], endo[:λ_f_t]]  =  1.
            Γ1[regime][eq[:eq_λ_f], endo[:λ_f_t]]  =  m[:ρ_λ_f]
            Γ1[regime][eq[:eq_λ_f], endo[:λ_f_t1]] = -m[:η_λ_f]
        Ψ[regime][eq[:eq_λ_f], exo[:λ_f_sh]]   =  1.

            Γ0[regime][eq[:eq_λ_f1], endo[:λ_f_t1]] = 1.
        Ψ[regime][eq[:eq_λ_f1], exo[:λ_f_sh]]   = 1.

    # Wage mark-up shock
            Γ0[regime][eq[:eq_λ_w], endo[:λ_w_t]]  =  1.
            Γ1[regime][eq[:eq_λ_w], endo[:λ_w_t]]  =  m[:ρ_λ_w]
            Γ1[regime][eq[:eq_λ_w], endo[:λ_w_t1]] = -m[:η_λ_w]
        Ψ[regime][eq[:eq_λ_w], exo[:λ_w_sh]]   =  1.

            Γ0[regime][eq[:eq_λ_w1], endo[:λ_w_t1]] = 1.
       Ψ[regime][eq[:eq_λ_w1], exo[:λ_w_sh]]   = 1.

    # Monetary policy shock
            Γ0[regime][eq[:eq_rm], endo[:rm_t]] = 1.
            Γ1[regime][eq[:eq_rm], endo[:rm_t]] = m[:ρ_rm]
        Ψ[regime][eq[:eq_rm], exo[:rm_sh]]  = 1.

    # Anticipated policy shocks
    if n_anticipated_shocks(m) > 0

        # This section adds the anticipated shocks. There is one state for all the
        # anticipated shocks that will hit in a given period (i.e. rm_tl2 holds those that
        # will hit in two periods), and the equations are set up so that rm_tl2 last period
        # will feed into rm_tl1 this period (and so on for other numbers), and last period's
        # rm_tl1 will feed into the rm_t process (and affect the Taylor Rule this period).

                Γ1[regime][eq[:eq_rm], endo[:rm_tl1]]   = 1.
                Γ0[regime][eq[:eq_rml1], endo[:rm_tl1]] = 1.
        Ψ[regime][eq[:eq_rml1], exo[:rm_shl1]]  = 1.

        if n_anticipated_shocks(m) > 1
            for i = 2:n_anticipated_shocks(m)
                        Γ1[regime][eq[Symbol("eq_rml$(i-1)")], endo[Symbol("rm_tl$i")]] = 1.
                        Γ0[regime][eq[Symbol("eq_rml$i")], endo[Symbol("rm_tl$i")]] = 1.
                Ψ[regime][eq[Symbol("eq_rml$i")], exo[Symbol("rm_shl$i")]] = 1.
            end
        end
    end



    ### RATIONAL EXPECTATIONS ERRORS ###

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


    ### E(rk)

    # Sticky prices and wages
            Γ0[regime][eq[:eq_Erk], endo[:rk_t]]  = 1.
            Γ1[regime][eq[:eq_Erk], endo[:Erk_t]] = 1.
        Π[regime][eq[:eq_Erk], ex[:Erk_sh]]   = 1.

    # Flexible prices and wages
            Γ0[regime][eq[:eq_Erk_f], endo[:rk_f_t]]  = 1.
            Γ1[regime][eq[:eq_Erk_f], endo[:Erk_f_t]] = 1.
        Π[regime][eq[:eq_Erk_f], ex[:Erk_f_sh]]   = 1.


    ### E(w)

    # Sticky prices and wages
            Γ0[regime][eq[:eq_Ew], endo[:w_t]]  = 1.
            Γ1[regime][eq[:eq_Ew], endo[:Ew_t]] = 1.
            Π[regime][eq[:eq_Ew], ex[:Ew_sh]]   = 1.
        end
        return (Γ0[1], Γ1[1], C[1], Ψ[1], Π[1]), (Γ0[2], Γ1[2], C[2], Ψ[2], Π[2])
end

   #= elseif subspec(m) =="ss51"
        for regime = 1:2
            endo = m.endogenous_states
            exo  = m.exogenous_shocks
            ex   = m.expected_shocks
            eq   = m.equilibrium_conditions

            Γ0[regime] = zeros(n_states(m), n_states(m))
            Γ1[regime] = zeros(n_states(m), n_states(m))
            C[regime]  = zeros(n_states(m))
            Ψ[regime]  = zeros(n_states(m), n_shocks_exogenous(m))
            Π[regime]  = zeros(n_states(m), n_shocks_expectational(m))

            if regime == 1
                S′′    = m[:S′′]
                σ_c    = m[:σ_c]
                h      = m[:h]
                λ_w    = m[:λ_w]
                ζ_w    = m[:ζ_w]
                ϵ_w    = m[:ϵ_w]
                ν_l    = m[:ν_l]
                δ      = m[:δ]
                ζ_p    = m[:ζ_p]
                ϵ_p    = m[:ϵ_p]
                ι_w    = m[:ι_w]
                ι_p    = m[:ι_p]
                ppsi   = m[:ppsi]
                Φ      = m[:Φ]
                ψ1     = m[:ψ1]
                ρ      = m[:ρ]
                ψ2     = m[:ψ2]
                ψ3     = m[:ψ3]
                π_star = m[:π_star]
                β      = m[:β]
                Lmean  = m[:Lmean]
                γ      = m[:γ]
                α      = m[:α]
                g_star = m[:g_star]
                ρ_z    = m[:ρ_z]
                ρ_b    = m[:ρ_b]
                ρ_g    = m[:ρ_g]
                ρ_μ    = m[:ρ_μ]
                ρ_rm   = m[:ρ_rm]
                ρ_λ_f  = m[:ρ_λ_f]
                ρ_λ_w  = m[:ρ_λ_w]
                η_λ_f  = m[:η_λ_f]
                η_λ_w  = m[:η_λ_w]
                η_gz   = m[:η_gz]
                σ_z    = m[:σ_z]
                σ_b    = m[:σ_b]
                σ_g    = m[:σ_g]
                σ_μ    = m[:σ_μ]
                σ_rm   = m[:σ_rm]
                σ_λ_f  = m[:σ_λ_f]
                σ_λ_w  = m[:σ_λ_w]
                Rstarn = m[:Rstarn]
                rkstar = m[:rkstar]
                wstar  = m[:wstar]
                i_k    = m[:i_k]
                l_k    = m[:l_k]
                k_y    = m[:k_y]
                i_y    = m[:i_y]
                c_y    = m[:c_y]
                u_y    = m[:u_y]
                wl_c   = m[:wl_c]
            elseif regime == 2
                S′′    = m[:S′′_r2]
                σ_c    = m[:σ_c_r2]
                h      = m[:h_r2]
                λ_w    = m[:λ_w_r2]
                ζ_w    = m[:ζ_w_r2]
                ϵ_w    = m[:ϵ_w_r2]
                ν_l    = m[:ν_l_r2]
                δ      = m[:δ_r2]
                ζ_p    = m[:ζ_p_r2]
                ϵ_p    = m[:ϵ_p_r2]
                ι_w    = m[:ι_w_r2]
                ι_p    = m[:ι_p_r2]
                ppsi   = m[:ppsi_r2]
                Φ      = m[:Φ_r2]
                ψ1     = m[:ψ1_r2]
                ρ      = m[:ρ_r2]
                ψ2     = m[:ψ2_r2]
                ψ3     = m[:ψ3_r2]
                π_star = m[:π_star_r2]
                β      = m[:β_r2]
                Lmean  = m[:Lmean_r2]
                γ      = m[:γ_r2]
                α      = m[:α_r2]
                g_star = m[:g_star_r2]
                ρ_z    = m[:ρ_z_r2]
                ρ_b    = m[:ρ_b_r2]
                ρ_g    = m[:ρ_g_r2]
                ρ_μ    = m[:ρ_μ_r2]
                ρ_rm   = m[:ρ_rm_r2]
                ρ_λ_f  = m[:ρ_λ_f_r2]
                ρ_λ_w  = m[:ρ_λ_w_r2]
                η_λ_f  = m[:η_λ_f_r2]
                η_λ_w  = m[:η_λ_w_r2]
                η_gz   = m[:η_gz_r2]
                σ_z    = m[:σ_z_r2]
                σ_b    = m[:σ_b_r2]
                σ_g    = m[:σ_g_r2]
                σ_μ    = m[:σ_μ_r2]
                σ_rm   = m[:σ_rm_r2]
                σ_λ_f  = m[:σ_λ_f_r2]
                σ_λ_w  = m[:σ_λ_w_r2]
                Rstarn = m[:Rstarn_r2]
                rkstar = m[:rkstar_r2]
                wstar  = m[:wstar_r2]
                i_k    = m[:i_k_r2]
                l_k    = m[:l_k_r2]
                k_y    = m[:k_y_r2]
                i_y    = m[:i_y_r2]
                c_y    = m[:c_y_r2]
                u_y    = m[:u_y_r2]
                wl_c   = m[:wl_c_r2]
            end
            ### ENDOGENOUS STATES ###

            ### 1. Consumption Euler Equation

            #TRYING TO MATCH TO FORTRAN
            #γ] = γ]/100 + 1
            #β] = 1/(β]/100+1)
            #π_star] = π_star]/100+1

            # Sticky prices and wages
            Γ0[regime][eq[:eq_euler], endo[:c_t]]    = 1.
            Γ0[regime][eq[:eq_euler], endo[:Ec_t]]   = -1/(1 + h/γ)
            Γ0[regime][eq[:eq_euler], endo[:L_t]]    = -(σ_c - 1)*wl_c/(σ_c*(1 + h/γ))
            Γ0[regime][eq[:eq_euler], endo[:EL_t]]   = (σ_c - 1)*wl_c/(σ_c*(1 + h/γ))
            Γ0[regime][eq[:eq_euler], endo[:R_t]]    = (1 - h/γ)/(σ_c*(1 + h/γ))
        Γ0[regime][eq[:eq_euler], endo[:Eπ_t]]   = -(1 - h/γ)/(σ_c*(1 + h/γ))
        Γ0[regime][eq[:eq_euler], endo[:b_t]]    = -1.
        Γ1[regime][eq[:eq_euler], endo[:c_t]]    = (h/γ)/(1 + h/γ)

        # Flexible prices and wages
        Γ0[regime][eq[:eq_euler_f], endo[:c_f_t]]  = 1.
        Γ0[regime][eq[:eq_euler_f], endo[:Ec_f_t]] = -1/(1 + h/γ)
        Γ0[regime][eq[:eq_euler_f], endo[:L_f_t]]  = -(σ_c - 1)*wl_c/(σ_c*(1 + h/γ))
        Γ0[regime][eq[:eq_euler_f], endo[:EL_f_t]] = (σ_c - 1)*wl_c/(σ_c*(1 + h/γ))
        Γ0[regime][eq[:eq_euler_f], endo[:r_f_t]]  = (1 - h/γ)/(σ_c*(1 + h/γ))
        Γ0[regime][eq[:eq_euler_f], endo[:b_t]]    = -1.
        Γ1[regime][eq[:eq_euler_f], endo[:c_f_t]]  = (h/γ)/(1 + h/γ)


        ### 2. Investment Euler Equation

        # Sticky prices and wages
        Γ0[regime][eq[:eq_inv], endo[:i_t]]  = 1.
        Γ0[regime][eq[:eq_inv], endo[:Ei_t]] = -β*γ^(1 - σ_c)/(1 + β*γ^(1 - σ_c))
        Γ0[regime][eq[:eq_inv], endo[:qk_t]] = -1 / (S′′ * γ^2 * (1 + β*γ^(1 - σ_c)))
        Γ0[regime][eq[:eq_inv], endo[:μ_t]]  = -1.
        Γ1[regime][eq[:eq_inv], endo[:i_t]]  = 1/(1 + β*γ^(1 - σ_c))

        # Flexible prices and wages
        Γ0[regime][eq[:eq_inv_f], endo[:i_f_t]]  = 1.
        Γ0[regime][eq[:eq_inv_f], endo[:Ei_f_t]] = -β*γ^(1 - σ_c)/(1 + β*γ^(1 - σ_c))
        Γ0[regime][eq[:eq_inv_f], endo[:qk_f_t]] = -1 / (S′′ * γ^2 * (1 + β*γ^(1 - σ_c)))
        Γ0[regime][eq[:eq_inv_f], endo[:μ_t]]    = -1.
        Γ1[regime][eq[:eq_inv_f], endo[:i_f_t]]  = 1/(1 + β*γ^(1 - σ_c))


        ### 3. Value of Capital

        # Sticky prices and wages

        Γ0[regime][eq[:eq_capval],endo[:qk_t]]  = 1.
        Γ0[regime][eq[:eq_capval],endo[:Eqk_t]] = -β * γ^(-σ_c) * (1 - δ)
        Γ0[regime][eq[:eq_capval],endo[:Erk_t]] = -(1 - β * γ^(-σ_c) * (1 - δ))
        Γ0[regime][eq[:eq_capval],endo[:R_t]]   = 1.
        Γ0[regime][eq[:eq_capval],endo[:Eπ_t]]  = -1.
        Γ0[regime][eq[:eq_capval],endo[:b_t]]   = -1 / ((1-h/γ) / ((1+h/γ)* σ_c))

        # Flexible prices and wages
        Γ0[regime][eq[:eq_capval_f], endo[:qk_f_t]]  = 1.
        Γ0[regime][eq[:eq_capval_f], endo[:Eqk_f_t]] = -β * γ^(-σ_c) * (1 - δ)
        Γ0[regime][eq[:eq_capval_f], endo[:Erk_f_t]] = -(1 - β * γ^(-σ_c) * (1 - δ))
        Γ0[regime][eq[:eq_capval_f], endo[:r_f_t]]   = 1.
        Γ0[regime][eq[:eq_capval_f], endo[:b_t]]     = -1 / ((1-h/γ) / ((1+h/γ)* σ_c))


        ### 4. Aggregate Production Function

        # Sticky prices and wages
        Γ0[regime][eq[:eq_output], endo[:y_t]] =  1.
        Γ0[regime][eq[:eq_output], endo[:k_t]] = -Φ*α
        Γ0[regime][eq[:eq_output], endo[:L_t]] = -Φ*(1 - α)
        Γ0[regime][eq[:eq_output], endo[:z_t]] = -Φ

        # Flexible prices and wages
        Γ0[regime][eq[:eq_output_f], endo[:y_f_t]] =  1.
        Γ0[regime][eq[:eq_output_f], endo[:k_f_t]] = -Φ*α
        Γ0[regime][eq[:eq_output_f], endo[:L_f_t]] = -Φ*(1 - α)
        Γ0[regime][eq[:eq_output_f], endo[:z_t]]   = -Φ


    ### 5. Capital Utilization

    # Sticky prices and wages
        Γ0[regime][eq[:eq_caputl], endo[:k_t]]    =  1.
        Γ1[regime][eq[:eq_caputl], endo[:kbar_t]] =  1.
        Γ0[regime][eq[:eq_caputl], endo[:u_t]]    = -1.

    # Flexible prices and wages
        Γ0[regime][eq[:eq_caputl_f], endo[:k_f_t]]    =  1.
        Γ1[regime][eq[:eq_caputl_f], endo[:kbar_f_t]] =  1.
        Γ0[regime][eq[:eq_caputl_f], endo[:u_f_t]]    = -1.



    ### 6. Rental Rate of Capital

    # Sticky prices and wages
        Γ0[regime][eq[:eq_capsrv], endo[:u_t]]  = 1.
        Γ0[regime][eq[:eq_capsrv], endo[:rk_t]] = -(1 - ppsi)/ppsi

    # Flexible prices and wages
        Γ0[regime][eq[:eq_capsrv_f], endo[:u_f_t]]  = 1.
        Γ0[regime][eq[:eq_capsrv_f], endo[:rk_f_t]] = -(1 - ppsi)/ppsi



    ### 7. Evolution of Capital

    # Sticky prices and wages
        Γ0[regime][eq[:eq_capev], endo[:kbar_t]] = 1.
        Γ0[regime][eq[:eq_capev], endo[:μ_t]]    = -(1 - (1 - δ)/γ) * γ^2 * S′′ * (1 + β*γ^(1 - σ_c))
        Γ0[regime][eq[:eq_capev], endo[:i_t]]    = -(1 - (1 - δ)/γ)
        Γ1[regime][eq[:eq_capev], endo[:kbar_t]] = (1 - δ)/γ

    # Flexible prices and wages
        Γ0[regime][eq[:eq_capev_f], endo[:kbar_f_t]] = 1.
        Γ0[regime][eq[:eq_capev_f], endo[:μ_t]]      = -(1 - (1 - δ)/γ) * γ^2 * S′′ * (1 + β*γ^(1 - σ_c))
        Γ0[regime][eq[:eq_capev_f], endo[:i_f_t]]    = -(1 - (1 - δ)/γ)
        Γ1[regime][eq[:eq_capev_f], endo[:kbar_f_t]] = (1 - δ)/γ


    ### 8. Price Markup

    # Sticky prices and wages
        Γ0[regime][eq[:eq_mkupp], endo[:mc_t]] =  1.
        Γ0[regime][eq[:eq_mkupp], endo[:k_t]]  = -α
        Γ0[regime][eq[:eq_mkupp], endo[:L_t]]  =  α
        Γ0[regime][eq[:eq_mkupp], endo[:z_t]]  = -1.
        Γ0[regime][eq[:eq_mkupp], endo[:w_t]]  =  1.

    # Flexible prices and wages
        Γ0[regime][eq[:eq_mkupp_f], endo[:k_f_t]] = -α
        Γ0[regime][eq[:eq_mkupp_f], endo[:L_f_t]] =  α
        Γ0[regime][eq[:eq_mkupp_f], endo[:z_t]]   = -1.
        Γ0[regime][eq[:eq_mkupp_f], endo[:w_f_t]] =  1.


    ### 9. Phillips Curve

    # Sticky prices and wages
        Γ0[regime][eq[:eq_phlps], endo[:π_t]]   = 1.
        Γ0[regime][eq[:eq_phlps], endo[:Eπ_t]]  = -β * γ^(1 - σ_c) / (1 + β*γ^(1 - σ_c)*ι_p)
        Γ0[regime][eq[:eq_phlps], endo[:mc_t]]  = ((1 - ζ_p)*(1 - ζ_p*β*γ^(1 - σ_c))) /
            (ζ_p*((Φ- 1)*ϵ_p + 1)*(1 + β*γ^(1 - σ_c)*ι_p))

        if regime == 1
            Γ0[regime][eq[:eq_phlps], endo[:λ_f_t]] = -1.
        elseif regime == 2
            κnum = ((1 - ζ_p*β*exp((1 - σ_c)))* #zstar]))*
            (1 - ζ_p))/(ζ_p*((Φ- 1)*ϵ_p + 1))/(1 + ι_p*β*exp((1 - σ_c))) #*zstar]))         # kappa numerator
                fix_ζ_p = ζ_p
                κden = ((1 - fix_ζ_p*β*exp((1 - σ_c)))* #zstar]))*
                        (1 - fix_ζ_p))/(fix_ζ_p*((Φ- 1)*ϵ_p + 1))/(1 + ι_p*β*exp((1 - σ_c))) #*zstar]))         # kappa denominator
                Γ0[regime][eq[:eq_phlps], endo[:λ_f_t]] = -κnum / κden
        else
            Γ0[regime][eq[:eq_phlps], endo[:λ_f_t]] = -1.
        end
        Γ1[regime][eq[:eq_phlps], endo[:π_t]]   = ι_p / (1 + β*γ^(1 - σ_c)*ι_p)

    # Flexible prices and wages not necessary

    ### 10. Rental Rate of Capital

    # Sticky prices and wages
        Γ0[regime][eq[:eq_caprnt], endo[:rk_t]] =  1.
        Γ0[regime][eq[:eq_caprnt], endo[:k_t]]  =  1.
        Γ0[regime][eq[:eq_caprnt], endo[:L_t]]  = -1.
        Γ0[regime][eq[:eq_caprnt], endo[:w_t]]  = -1.

    # Flexible prices and wages
        Γ0[regime][eq[:eq_caprnt_f], endo[:rk_f_t]] =  1.
        Γ0[regime][eq[:eq_caprnt_f], endo[:k_f_t]]  =  1.
        Γ0[regime][eq[:eq_caprnt_f], endo[:L_f_t]]  = -1.
        Γ0[regime][eq[:eq_caprnt_f], endo[:w_f_t]]  = -1.


    ### 11. Wage Markup (Marginal Substitution)

    # Sticky prices and wages
        Γ0[regime][eq[:eq_msub], endo[:μ_ω_t]] =  1.
        Γ0[regime][eq[:eq_msub], endo[:w_t]]   = -1.
        Γ0[regime][eq[:eq_msub], endo[:L_t]]   = ν_l
        Γ0[regime][eq[:eq_msub], endo[:c_t]]   = 1/(1 - h/γ)
        Γ1[regime][eq[:eq_msub], endo[:c_t]]   = (h/γ) / (1 - h/γ)

    # Flexible prices and wages
        Γ0[regime][eq[:eq_msub_f], endo[:w_f_t]] = -1.
        Γ0[regime][eq[:eq_msub_f], endo[:L_f_t]] = ν_l
        Γ0[regime][eq[:eq_msub_f], endo[:c_f_t]] = 1/(1 - h/γ)
        Γ1[regime][eq[:eq_msub_f], endo[:c_f_t]] = (h/γ) / (1 - h/γ)


    ### 12. Evolution of Wages

    # Sticky prices and wages
        Γ0[regime][eq[:eq_wage], endo[:w_t]]    = 1
        Γ0[regime][eq[:eq_wage], endo[:Ew_t]]   = -(1 - 1/(1 + β*γ^(1 - σ_c)))
        Γ0[regime][eq[:eq_wage], endo[:Eπ_t]]   = -(1 - 1/(1 + β*γ^(1 - σ_c)))
        Γ0[regime][eq[:eq_wage], endo[:π_t]]    = (1 + β*γ^(1 - σ_c)*ι_w) / (1 + β*γ^(1 - σ_c))
        Γ0[regime][eq[:eq_wage], endo[:μ_ω_t]]  = ((1 - β*γ^(1 - σ_c)*ζ_w) * (1 - ζ_w)) /
        ((1 + β*γ^(1 - σ_c)) * (ζ_w*((λ_w - 1)*ϵ_w + 1)))
        Γ0[regime][eq[:eq_wage], endo[:λ_w_t]]  = -1.
        Γ1[regime][eq[:eq_wage], endo[:w_t]]    = 1/(1 + β*γ^(1 - σ_c))
        Γ1[regime][eq[:eq_wage], endo[:π_t]]    = ι_w/(1 + β*γ^(1 - σ_c))

    # Flexible prices and wages not necessary


    ### 13. Monetary Policy Rule

    # Sticky prices and wages
        Γ0[regime][eq[:eq_mp], endo[:R_t]]   = 1.
        Γ1[regime][eq[:eq_mp], endo[:R_t]]   = ρ
            Γ0[regime][eq[:eq_mp], endo[:π_t]]   = -(1 - ρ)*ψ1
            Γ0[regime][eq[:eq_mp], endo[:y_t]]   = -(1 - ρ)*ψ2 - ψ3
            Γ1[regime][eq[:eq_mp], endo[:y_t]]   = -ψ3
            Γ0[regime][eq[:eq_mp], endo[:y_f_t]] = (1 - ρ)*ψ2 + ψ3
            Γ1[regime][eq[:eq_mp], endo[:y_f_t]] = ψ3
            Γ0[regime][eq[:eq_mp], endo[:rm_t]]  = -1.

    # Flexible prices and wages not necessary


    ### 14. Resource Constraint

    # Sticky prices and wages
        Γ0[regime][eq[:eq_res], endo[:y_t]] =  1.
        Γ0[regime][eq[:eq_res], endo[:c_t]] = -c_y
        Γ0[regime][eq[:eq_res], endo[:i_t]] = -i_y
        Γ0[regime][eq[:eq_res], endo[:u_t]] = -u_y
        Γ0[regime][eq[:eq_res], endo[:g_t]] = -1.

    # Flexible prices and wages
        Γ0[regime][eq[:eq_res_f], endo[:y_f_t]] = 1.
        Γ0[regime][eq[:eq_res_f], endo[:c_f_t]] = -c_y
        Γ0[regime][eq[:eq_res_f], endo[:i_f_t]] = -i_y
        Γ0[regime][eq[:eq_res_f], endo[:u_f_t]] = -u_y
        Γ0[regime][eq[:eq_res_f], endo[:g_t]]   = -1.


    ### EXOGENOUS PROCESSES ###

    # Government spending
        Γ0[regime][eq[:eq_g], endo[:g_t]] = 1.
        Γ1[regime][eq[:eq_g], endo[:g_t]] = ρ_g
        Ψ[regime][eq[:eq_g], exo[:g_sh]]  = 1.
        Ψ[regime][eq[:eq_g], exo[:z_sh]]  = η_gz

    # Asset shock
        Γ0[regime][eq[:eq_b], endo[:b_t]] = 1.
        Γ1[regime][eq[:eq_b], endo[:b_t]] = ρ_b
        Ψ[regime][eq[:eq_b], exo[:b_sh]]  = 1.

    # Investment-specific technology
        Γ0[regime][eq[:eq_μ], endo[:μ_t]] = 1.
        Γ1[regime][eq[:eq_μ], endo[:μ_t]] = ρ_μ
        Ψ[regime][eq[:eq_μ], exo[:μ_sh]]  = 1.

    # Neutral technology
        Γ0[regime][eq[:eq_z], endo[:z_t]] = 1.
        Γ1[regime][eq[:eq_z], endo[:z_t]] = ρ_z
        Ψ[regime][eq[:eq_z], exo[:z_sh]]  = 1.

    # Price mark-up shock
        Γ0[regime][eq[:eq_λ_f], endo[:λ_f_t]]  =  1.
        Γ1[regime][eq[:eq_λ_f], endo[:λ_f_t]]  =  ρ_λ_f
        Γ1[regime][eq[:eq_λ_f], endo[:λ_f_t1]] = -η_λ_f
        Ψ[regime][eq[:eq_λ_f], exo[:λ_f_sh]]   =  1.

        Γ0[regime][eq[:eq_λ_f1], endo[:λ_f_t1]] = 1.
        Ψ[regime][eq[:eq_λ_f1], exo[:λ_f_sh]]   = 1.

    # Wage mark-up shock
        Γ0[regime][eq[:eq_λ_w], endo[:λ_w_t]]  =  1.
        Γ1[regime][eq[:eq_λ_w], endo[:λ_w_t]]  =  ρ_λ_w
        Γ1[regime][eq[:eq_λ_w], endo[:λ_w_t1]] = -η_λ_w
        Ψ[regime][eq[:eq_λ_w], exo[:λ_w_sh]]   =  1.

        Γ0[regime][eq[:eq_λ_w1], endo[:λ_w_t1]] = 1.
       Ψ[regime][eq[:eq_λ_w1], exo[:λ_w_sh]]   = 1.

    # Monetary policy shock
        Γ0[regime][eq[:eq_rm], endo[:rm_t]] = 1.
        Γ1[regime][eq[:eq_rm], endo[:rm_t]] = ρ_rm
        Ψ[regime][eq[:eq_rm], exo[:rm_sh]]  = 1.

    # Anticipated policy shocks
    if n_anticipated_shocks(m) > 0

        # This section adds the anticipated shocks. There is one state for all the
        # anticipated shocks that will hit in a given period (i.e. rm_tl2 holds those that
        # will hit in two periods), and the equations are set up so that rm_tl2 last period
        # will feed into rm_tl1 this period (and so on for other numbers), and last period's
        # rm_tl1 will feed into the rm_t process (and affect the Taylor Rule this period).

            Γ1[regime][eq[:eq_rm], endo[:rm_tl1]]   = 1.
            Γ0[regime][eq[:eq_rml1], endo[:rm_tl1]] = 1.
        Ψ[regime][eq[:eq_rml1], exo[:rm_shl1]]  = 1.

        if n_anticipated_shocks(m) > 1
            for i = 2:n_anticipated_shocks(m)
                    Γ1[regime][eq[Symbol("eq_rml$(i-1)")], endo[Symbol("rm_tl$i")]] = 1.
                    Γ0[regime][eq[Symbol("eq_rml$i")], endo[Symbol("rm_tl$i")]] = 1.
                Ψ[regime][eq[Symbol("eq_rml$i")], exo[Symbol("rm_shl$i")]] = 1.
            end
        end
    end



    ### RATIONAL EXPECTATIONS ERRORS ###

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


    ### E(rk)

    # Sticky prices and wages
        Γ0[regime][eq[:eq_Erk], endo[:rk_t]]  = 1.
        Γ1[regime][eq[:eq_Erk], endo[:Erk_t]] = 1.
        Π[regime][eq[:eq_Erk], ex[:Erk_sh]]   = 1.

    # Flexible prices and wages
        Γ0[regime][eq[:eq_Erk_f], endo[:rk_f_t]]  = 1.
        Γ1[regime][eq[:eq_Erk_f], endo[:Erk_f_t]] = 1.
        Π[regime][eq[:eq_Erk_f], ex[:Erk_f_sh]]   = 1.


    ### E(w)

    # Sticky prices and wages
        Γ0[regime][eq[:eq_Ew], endo[:w_t]]  = 1.
        Γ1[regime][eq[:eq_Ew], endo[:Ew_t]] = 1.
        Π[regime][eq[:eq_Ew], ex[:Ew_sh]]   = 1.
    end
    end
    return (Γ0[1], Γ1[1], C[1], Ψ[1], Π[1]), (Γ0[2], Γ1[2], C[2], Ψ[2], Π[2])
end =#
