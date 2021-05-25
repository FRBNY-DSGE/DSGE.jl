"""
```
pseudo_measurement(m::Model1002{T}, TTT::Matrix{T}, RRR::Matrix{T},
    CCC::Vector{T}; reg::Int = 1) where {T<:AbstractFloat}
```

Assign pseudo-measurement equation (a linear combination of states):

```
x_t = ZZ_pseudo*s_t + DD_pseudo
```
"""
function pseudo_measurement(m::Model1002{T},
                            TTT::Matrix{T},
                            RRR::Matrix{T},
                            CCC::Vector{T};
                            reg::Int = 1,
                            TTTs::Vector{<: AbstractMatrix{T}} = Matrix{T}[],
                            CCCs::Vector{<: AbstractVector{T}} = Vector{T}[],
                            information_set::UnitRange = reg:reg,
                            memo::Union{ForwardMultipleExpectationsMemo, Nothing} = nothing) where {T <: AbstractFloat}

    endo      = m.endogenous_states
    endo_addl = m.endogenous_states_augmented
    pseudo    = m.pseudo_observables

    # Initialize pseudo ZZ and DD matrices
    _n_states = n_states_augmented(m)
    _n_pseudo = n_pseudo_observables(m)

    ZZ_pseudo = zeros(_n_pseudo, _n_states)
    DD_pseudo = zeros(_n_pseudo)

    # Set parameters
    for para in m.parameters
        if !isempty(para.regimes)
            if (haskey(get_settings(m), :model2para_regime) ? haskey(get_setting(m, :model2para_regime), para.key) : false)
                ModelConstructors.toggle_regime!(para, reg, get_setting(m, :model2para_regime)[para.key])
            else
                ModelConstructors.toggle_regime!(para, reg)
            end
        end
    end

    # Set up for calculating k-periods ahead expectations and expected sums
    permanent_t = length(information_set[findfirst(x -> x == reg, information_set):end]) - 1 + reg
    if information_set[1] == information_set[end]
        # In this case, we do not need to pass the TTTs, CCCs in, so we redefine them as empty.
        # This step is also necessary to ensure the memo is properly used.
        TTTs = Matrix{T}[]
        CCCs = Vector{T}[]

        memo = nothing # see measurement
        # TODO: maybe instead of emptying TTTs, CCCs, we add
        # a step to recompute the memo since we will still likely be recalculating
        # products/powers of TTT multiple times that could be pre-computed
    end
    use_fwd_exp_sum = haskey(get_settings(m), :use_forward_expected_sum_memo) && get_setting(m, :use_forward_expected_sum_memo)
    use_fwd_exp     = haskey(get_settings(m), :use_forward_expectations_memo) && get_setting(m, :use_forward_expectations_memo)

    # Handle integrated series
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

    # Compute TTT^10, used for Expected10YearRateGap, Expected10YearRate, and Expected10YearNaturalRate
    TTT10, CCC10 = k_periods_ahead_expected_sums(TTT, CCC, TTTs, CCCs, reg, 40, permanent_t;
                                                 integ_series = integ_series,
                                                 memo = use_fwd_exp_sum ? memo : nothing)
    TTT10        = TTT10./ 40. # divide by 40 to average across 10 years
    CCC10        = CCC10 ./ 40.

    if get_setting(m, :add_laborproductivity_measurement)
        # Construct pseudo-obs from integrated states first
        ZZ_pseudo[pseudo[:laborproductivity], endo[:y_t]] = 1.
        ZZ_pseudo[pseudo[:laborproductivity], endo[:L_t]] = -1.
        # ZZ_pseudo[pseudo[:laborproductivity], endo_addl[:cum_z_t]] = 1.
        DD_pseudo[pseudo[:laborproductivity]] = 100. * log(m[:ystar] / m[:Lstar])

        # Remove integrated states (e.g. states w/unit roots)
        # RRR and CCC aren't used, so we don't do anything with them
    end

    if get_setting(m, :add_nominalgdp_level)
        ZZ_pseudo[pseudo[:NominalGDPLevel], endo_addl[:cum_y_t]]     = 1.
        ZZ_pseudo[pseudo[:NominalGDPLevel], endo_addl[:cum_z_t]]     = 1.
        ZZ_pseudo[pseudo[:NominalGDPLevel], endo_addl[:cum_e_gdp_t]] = 1.
        ZZ_pseudo[pseudo[:NominalGDPLevel], endo_addl[:cum_π_t]]     = 1.
    end

    if get_setting(m, :add_nominalgdp_growth)
        ZZ_pseudo[pseudo[:NominalGDPGrowth], endo[:y_t]]           = 1.
        ZZ_pseudo[pseudo[:NominalGDPGrowth], endo_addl[:y_t1]]     = -1.
        ZZ_pseudo[pseudo[:NominalGDPGrowth], endo_addl[:z_t]]      = 1.
        ZZ_pseudo[pseudo[:NominalGDPGrowth], endo_addl[:e_gdp_t]]  = 1.
        ZZ_pseudo[pseudo[:NominalGDPGrowth], endo_addl[:e_gdp_t1]] = -m[:me_level]
        ZZ_pseudo[pseudo[:NominalGDPGrowth], endo[:π_t]]           = 1.
        DD_pseudo[pseudo[:NominalGDPGrowth]]                       = 100. * (exp(m[:z_star] - 1.) + (m[:π_star] - 1.))
    end

    if get_setting(m, :add_cumulative)
        ZZ_pseudo[pseudo[:AccumOutputGap], endo_addl[:cum_y_t]]   = 1.
        ZZ_pseudo[pseudo[:AccumOutputGap], endo_addl[:cum_y_f_t]] = -1.

        ZZ_pseudo[pseudo[:GDPLevel], endo_addl[:cum_y_t]]     = 1.
        ZZ_pseudo[pseudo[:GDPLevel], endo_addl[:cum_z_t]]     = 1.
        ZZ_pseudo[pseudo[:GDPLevel], endo_addl[:cum_e_gdp_t]] = 1.

        ZZ_pseudo[pseudo[:FlexibleGDPLevel], endo_addl[:cum_y_f_t]] = 1.
        ZZ_pseudo[pseudo[:FlexibleGDPLevel], endo_addl[:cum_z_t]]   = 1.

        ZZ_pseudo[pseudo[:ConsumptionLevel], endo_addl[:cum_c_t]]     = 1.
        ZZ_pseudo[pseudo[:ConsumptionLevel], endo_addl[:cum_z_t]]     = 1.

        ZZ_pseudo[pseudo[:FlexibleConsumptionLevel], endo_addl[:cum_c_f_t]] = 1.
        ZZ_pseudo[pseudo[:FlexibleConsumptionLevel], endo_addl[:cum_z_t]]   = 1.

        ZZ_pseudo[pseudo[:InvestmentLevel], endo_addl[:cum_i_t]]     = 1.
        ZZ_pseudo[pseudo[:InvestmentLevel], endo_addl[:cum_z_t]]     = 1.

        ZZ_pseudo[pseudo[:FlexibleInvestmentLevel], endo_addl[:cum_i_f_t]] = 1.
        ZZ_pseudo[pseudo[:FlexibleInvestmentLevel], endo_addl[:cum_z_t]]   = 1.
    end

    if get_setting(m, :add_flexible_price_growth)
        ZZ_pseudo[pseudo[:FlexibleGDPGrowth], endo[:y_f_t]]       = 1.
        ZZ_pseudo[pseudo[:FlexibleGDPGrowth], endo_addl[:y_f_t1]] = -1.
        ZZ_pseudo[pseudo[:FlexibleGDPGrowth], endo[:z_t]]         = 1.
        DD_pseudo[pseudo[:FlexibleGDPGrowth]]                     = 100. * (exp(m[:z_star]) - 1.)

        ZZ_pseudo[pseudo[:FlexibleInvestmentGrowth], endo[:i_f_t]]       = 1.
        ZZ_pseudo[pseudo[:FlexibleInvestmentGrowth], endo_addl[:i_f_t1]] = -1.
        ZZ_pseudo[pseudo[:FlexibleInvestmentGrowth], endo[:z_t]]         = 1.
        DD_pseudo[pseudo[:FlexibleInvestmentGrowth]]                     = 100. * (exp(m[:z_star]) - 1.)

        ZZ_pseudo[pseudo[:FlexibleConsumptionGrowth], endo[:c_f_t]]       = 1.
        ZZ_pseudo[pseudo[:FlexibleConsumptionGrowth], endo_addl[:c_f_t1]] = -1.
        ZZ_pseudo[pseudo[:FlexibleConsumptionGrowth], endo[:z_t]]         = 1.
        DD_pseudo[pseudo[:FlexibleConsumptionGrowth]]                     = 100. * (exp(m[:z_star]) - 1.)
    end

    ##########################################################
    ## PSEUDO-OBSERVABLE EQUATIONS
    ##########################################################

    ## Output
    ZZ_pseudo[pseudo[:y_t],endo[:y_t]] = 1.

    ## Flexible Output
    ZZ_pseudo[pseudo[:y_f_t],endo[:y_f_t]] = 1.

    ## Pseudo GDP Growth
    if haskey(m.settings, :add_pseudo_gdp)
        if get_setting(m, :add_pseudo_gdp) && subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63", "ss64", "ss65", "ss66", "ss67", "ss68", "ss69", "ss70", "ss71", "ss72", "ss73", "ss74", "ss75", "ss76", "ss77", "ss78", "ss79", "ss80", "ss81", "ss82", "ss83"]
            ZZ_pseudo[pseudo[:PseudoGDP], endo[:y_t]]          = 1.0
            ZZ_pseudo[pseudo[:PseudoGDP], endo_addl[:y_t1]]     = -1.0
            ZZ_pseudo[pseudo[:PseudoGDP], endo[:z_t]]          = 1.0
            ZZ_pseudo[pseudo[:PseudoGDP], endo_addl[:e_gdp_t]]  = 1.0
            ZZ_pseudo[pseudo[:PseudoGDP], endo_addl[:e_gdp_t1]] = -m[:me_level]
            DD_pseudo[pseudo[:PseudoGDP]]                      = 100*(exp(m[:z_star])-1)
        end
    end

    ## Pseudo Core PCE # TODO
    if haskey(m.settings, :add_pseudo_corepce)
        if get_setting(m, :add_pseudo_corepce) && subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63", "ss64", "ss65", "ss66", "ss67", "ss68", "ss69", "ss70", "ss71", "ss72", "ss73", "ss74", "ss75", "ss76", "ss77", "ss78", "ss79", "ss80", "ss81", "ss82", "ss83"]
            ZZ_pseudo[pseudo[:PseudoCorePCE], endo[:π_t]]              = 1.0
            ZZ_pseudo[pseudo[:PseudoCorePCE], endo_addl[:e_corepce_t]] = 1.0
            DD_pseudo[pseudo[:PseudoCorePCE]]                          = 100. * (m[:π_star] - 1.)
        end
    end

    ## Natural Rate
    ZZ_pseudo[pseudo[:NaturalRate],endo[:r_f_t]] = 1.
    DD_pseudo[pseudo[:NaturalRate]]              = 100.0*(m[:rstar]-1.0)

    ## π_t
    ZZ_pseudo[pseudo[:π_t],endo[:π_t]] = 1.
    DD_pseudo[pseudo[:π_t]]            = 100*(m[:π_star]-1);

    ## Output Gap
    ZZ_pseudo[pseudo[:OutputGap],endo[:y_t]] = 1;
    ZZ_pseudo[pseudo[:OutputGap],endo[:y_f_t]] = -1;

    ## Ex Ante Real Rate
    ZZ_pseudo[pseudo[:ExAnteRealRate],endo[:R_t]]  = 1;
    ZZ_pseudo[pseudo[:ExAnteRealRate],endo[:Eπ_t]] = -1;
    DD_pseudo[pseudo[:ExAnteRealRate]]             = m[:Rstarn] - 100*(m[:π_star]-1);

    ## Long Run Inflation
    ZZ_pseudo[pseudo[:LongRunInflation],endo[:π_star_t]] = 1.
    DD_pseudo[pseudo[:LongRunInflation]]                 = 100. *(m[:π_star]-1.)

    ## Marginal Cost
    ZZ_pseudo[pseudo[:MarginalCost],endo[:mc_t]] = 1.

    ## Wages
    ZZ_pseudo[pseudo[:Wages],endo[:w_t]] = 1.

    ## Flexible Wages
    ZZ_pseudo[pseudo[:FlexibleWages],endo[:w_f_t]] = 1.

    # ## b Wages
    # ZZ_pseudo[pseudo[:b_t],endo[:b_t]] = 1.

    ## Hours
    ZZ_pseudo[pseudo[:Hours],endo[:L_t]] = 1.

    ## Flexible Hours
    ZZ_pseudo[pseudo[:FlexibleHours],endo[:L_f_t]] = 1.

    ## z_t
    ZZ_pseudo[pseudo[:z_t], endo[:z_t]] = 1.

    ## Expected 10-Year Rate Gap
    # ZZ_pseudo[pseudo[:Expected10YearRateGap], :] = TTT10[endo[:R_t], :] - TTT10[endo[:r_f_t], :] - TTT10[endo[:Eπ_t], :]
    ZZ_pseudo[pseudo[:Expected10YearRateGap], :] = view(TTT10, endo[:R_t], :) - view(TTT10, endo[:r_f_t], :) - view(TTT10, endo[:Eπ_t], :)
    DD_pseudo[pseudo[:Expected10YearRateGap]]    = CCC10[endo[:R_t]] - CCC10[endo[:r_f_t]] - CCC10[endo[:Eπ_t]]

    ## Nominal FFR
    ZZ_pseudo[pseudo[:NominalFFR], endo[:R_t]] = 1.
    DD_pseudo[pseudo[:NominalFFR]] = m[:Rstarn]

    ## Expected 10-Year Interest Rate
    # ZZ_pseudo[pseudo[:Expected10YearRate], :] = TTT10[endo[:R_t], :]
    ZZ_pseudo[pseudo[:Expected10YearRate], :] = view(TTT10, endo[:R_t], :)
    DD_pseudo[pseudo[:Expected10YearRate]]    = m[:Rstarn] + CCC10[endo[:R_t]]

    ## Expected 10-Year Natural Rate
    # ZZ_pseudo[pseudo[:Expected10YearNaturalRate], :] = TTT10[endo[:r_f_t], :] + TTT10[endo[:Eπ_t], :]
    ZZ_pseudo[pseudo[:Expected10YearNaturalRate], :] = view(TTT10, endo[:r_f_t], :) + view(TTT10, endo[:Eπ_t], :)
    DD_pseudo[pseudo[:Expected10YearNaturalRate]]    = m[:Rstarn] + CCC10[endo[:r_f_t]] + CCC10[endo[:Eπ_t]]

    ## Expected Nominal Natural Rate
    ZZ_pseudo[pseudo[:ExpectedNominalNaturalRate], endo[:r_f_t]] = 1.
    ZZ_pseudo[pseudo[:ExpectedNominalNaturalRate], endo[:Eπ_t]]  = 1.
    DD_pseudo[pseudo[:ExpectedNominalNaturalRate]]               = m[:Rstarn]

    ## Nominal Rate Gap
    ZZ_pseudo[pseudo[:NominalRateGap], endo[:R_t]]   = 1.
    ZZ_pseudo[pseudo[:NominalRateGap], endo[:r_f_t]] = -1.
    ZZ_pseudo[pseudo[:NominalRateGap], endo[:Eπ_t]]  = -1.

    ## Labor Productivity Growth
    ZZ_pseudo[pseudo[:LaborProductivityGrowth], endo[:y_t]]           = 1.
    ZZ_pseudo[pseudo[:LaborProductivityGrowth], endo_addl[:y_t1]]     = -1.
    ZZ_pseudo[pseudo[:LaborProductivityGrowth], endo[:z_t]]           = 1.
    ZZ_pseudo[pseudo[:LaborProductivityGrowth], endo_addl[:e_gdp_t]]  = 1.
    ZZ_pseudo[pseudo[:LaborProductivityGrowth], endo_addl[:e_gdp_t1]] = -m[:me_level]
    ZZ_pseudo[pseudo[:LaborProductivityGrowth], endo[:L_t]]           = -1
    ZZ_pseudo[pseudo[:LaborProductivityGrowth], endo_addl[:L_t1]]     = 1.
    DD_pseudo[pseudo[:LaborProductivityGrowth]]                       = 100*(exp(m[:z_star]) - 1)

    ## u_t
    ZZ_pseudo[pseudo[:u_t], endo[:u_t]] = 1.

    ## Nominal Wage Growth
    if haskey(m.settings, :add_NominalWageGrowth)
        if get_setting(m, :add_NominalWageGrowth)
            ZZ_pseudo[pseudo[:NominalWageGrowth],endo[:w_t]] = 1.
            ZZ_pseudo[pseudo[:NominalWageGrowth],endo_addl[:w_t1]] = -1.
            ZZ_pseudo[pseudo[:NominalWageGrowth],endo[:z_t]] = 1.
            ZZ_pseudo[pseudo[:NominalWageGrowth],endo[:π_t]] = 1.
            DD_pseudo[pseudo[:NominalWageGrowth]]            = 100*(m[:π_star]-1) + 100*(exp(m[:z_star])-1)
        end
    end

    # ## i_f_t
    # ZZ_pseudo[pseudo[:i_f_t], endo[:i_f_t]] = 1.

    # ## R_t
    # ZZ_pseudo[pseudo[:R_t], endo[:R_t]] = 1.
    # DD_pseudo[pseudo[:R_t]] = 100.0*(m[:rstar]-1.0)

    # ## c_f_t
    # ZZ_pseudo[pseudo[:c_f_t], endo[:c_f_t]] = 1.

    # ## c_t
    # ZZ_pseudo[pseudo[:c_t], endo[:c_t]] = 1.

    # ## qk_f_t
    # ZZ_pseudo[pseudo[:qk_f_t], endo[:qk_f_t]] = 1.

    # ## k_f_t
    # ZZ_pseudo[pseudo[:k_f_t], endo[:k_f_t]] = 1.

    # ## r_f_t
    # ZZ_pseudo[pseudo[:r_f_t], endo[:r_f_t]] = 1.

    # ## kbar_f_t
    # ZZ_pseudo[pseudo[:kbar_f_t], endo[:kbar_f_t]] = 1.

    # ## u_f_t
    # ZZ_pseudo[pseudo[:u_f_t], endo[:u_f_t]] = 1.

    # ## rk_f_t
    # ZZ_pseudo[pseudo[:rk_f_t], endo[:rk_f_t]] = 1.

    # ## w_f_t
    # ZZ_pseudo[pseudo[:w_f_t], endo[:w_f_t]] = 1.

    # ## L_f_t
    # ZZ_pseudo[pseudo[:L_f_t], endo[:L_f_t]] = 1.

    # ## rktil_f_t
    # ZZ_pseudo[pseudo[:rktil_f_t], endo[:rktil_f_t]] = 1.

    # ## n_f_t
    # ZZ_pseudo[pseudo[:n_f_t], endo[:n_f_t]] = 1.

    ## labor share
    if haskey(m.settings, :add_laborshare_measurement)
        if get_setting(m, :add_laborshare_measurement)
            ZZ_pseudo[pseudo[:laborshare_t], endo[:w_t]] = 1.
            ZZ_pseudo[pseudo[:laborshare_t], endo[:L_t]] = 1.
            ZZ_pseudo[pseudo[:laborshare_t], endo[:y_t]] = -1.
            DD_pseudo[pseudo[:laborshare_t]] = 100. * log(m[:wstar] * m[:Lstar] / m[:ystar])
        end
    end

    ## Fundamental inflation related pseudo-obs
    if subspec(m) in ["ss13", "ss14", "ss15", "ss16", "ss17", "ss18", "ss19"]
        # Compute coefficient on Sinf
        betabar = exp((1-m[:σ_c] ) * m[:z_star]) * m[:β]
        κ = ((1 - m[:ζ_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))*
             (1 - m[:ζ_p]))/(m[:ζ_p]*((m[:Φ]- 1)*m[:ϵ_p] + 1))/
        (1 + m[:ι_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        κcoef = κ * (1 + m[:ι_p] * betabar)

        ZZ_pseudo[pseudo[:Sinf_t], endo_addl[:Sinf_t]] = 1.
        ZZ_pseudo[pseudo[:Sinf_w_coef_t], endo_addl[:Sinf_t]] = κcoef
        DD_pseudo[pseudo[:ι_p]] = m[:ι_p]
        ZZ_pseudo[pseudo[:πtil_t], endo_addl[:πtil_t]] = 1.
        DD_pseudo[pseudo[:πtil_t]] = 100 * (m[:π_star] - 1)
        ZZ_pseudo[pseudo[:e_tfp_t], endo_addl[:e_tfp_t]] = 1.
        if subspec(m) in ["ss14", "ss15", "ss16", "ss18", "ss19"]
            ZZ_pseudo[pseudo[:e_tfp_t1], endo_addl[:e_tfp_t1]] = 1.
        end
    end

    ## Exogenous processes
    if subspec(m) == "ss12"
        to_add = [:g_t, :b_t, :μ_t, :z_t, :λ_f_t, :λ_w_t, :rm_t, :σ_ω_t, :μ_e_t,
                  :γ_t, :π_star_t]
        to_add_addl = [:e_lr_t, :e_tfp_t, :e_gdpdef_t, :e_corepce_t, :e_gdp_t, :e_gdi_t]
        for i in to_add
            ZZ_pseudo[pseudo[i], endo[i]] = 1.
        end
        for i in to_add_addl
            ZZ_pseudo[pseudo[i], endo_addl[i]] = 1.
        end
    end

    if haskey(m.settings, :add_covid_pseudoobs)
        if get_setting(m, :add_covid_pseudoobs) && subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63", "ss64", "ss65", "ss66", "ss67", "ss68", "ss69", "ss70", "ss71", "ss72", "ss73", "ss74", "ss75", "ss76", "ss77", "ss78", "ss79", "ss80", "ss81", "ss82", "ss83"]
            ZZ_pseudo[pseudo[:ziid], endo[:ziid_t]] = 1.
            ZZ_pseudo[pseudo[:varphiiid], endo[:φ_t]] = 1.
            ZZ_pseudo[pseudo[:biidc], endo[:biidc_t]] = 1.
        end
    end

    if haskey(m.settings, :add_ztil)
        if get_setting(m, :add_ztil)
            ZZ_pseudo[pseudo[:ztil], endo[:ztil_t]] = 1.
        end
    end

    if haskey(m.settings, :add_zp)
        if get_setting(m, :add_zp)
            ZZ_pseudo[pseudo[:zp], endo[:zp_t]] = 1.
        end
    end

    if haskey(m.settings, :add_pgap)
        if get_setting(m, :add_pgap)
            ZZ_pseudo[pseudo[:pgap], endo[:pgap_t]] = 1.
        end
    end

    if haskey(m.settings, :add_ygap)
        if get_setting(m, :add_ygap)
            ZZ_pseudo[pseudo[:ygap], endo[:ygap_t]] = 1.
        end
    end

    if haskey(m.settings, :add_altpolicy_pgap)
        if get_setting(m, :add_altpolicy_pgap)
            ZZ_pseudo[pseudo[:pgap], endo[:pgap_t]] = 1.
        end
    end

    if haskey(m.settings, :add_altpolicy_ygap)
        if get_setting(m, :add_altpolicy_ygap)
            ZZ_pseudo[pseudo[:ygap], endo[:ygap_t]] = 1.
        end
    end

    if haskey(m.settings, :add_urhat)
        if get_setting(m, :add_urhat)
            ZZ_pseudo[pseudo[:urhat], endo[:L_t]] = -0.7366220966444673
            DD_pseudo[pseudo[:urhat]] = -32.257825316004364 + -0.7366220966444673 * m[:Lmean]
        end
    end

    if haskey(m.settings, :add_rw)
        if get_setting(m, :add_rw)
            ZZ_pseudo[pseudo[:rw], endo[:rw_t]]     = 1.
            ZZ_pseudo[pseudo[:Rref], endo[:Rref_t]] = 1.
            DD_pseudo[pseudo[:Rref]]                = m[:Rstarn]
        end
    end

    for para in m.parameters
        if !isempty(para.regimes)
            ModelConstructors.toggle_regime!(para, 1)
        end
    end

    return PseudoMeasurement(ZZ_pseudo, DD_pseudo)
end
