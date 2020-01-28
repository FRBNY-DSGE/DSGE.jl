"""
```
pseudo_measurement(m::Model1002{T},
    TTT::Matrix{T}, RRR::Matrix{T}, CCC::Vector{T}) where {T<:AbstractFloat}
```

Assign pseudo-measurement equation (a linear combination of states):

```
x_t = ZZ_pseudo*s_t + DD_pseudo
```
"""
function pseudo_measurement(m::Model1002{T},
                            TTT::Matrix{T},
                            RRR::Matrix{T},
                            CCC::Vector{T}) where {T<:AbstractFloat}

    endo      = m.endogenous_states
    endo_addl = m.endogenous_states_augmented
    pseudo    = m.pseudo_observables

    _n_states = n_states_augmented(m)
    _n_pseudo = n_pseudo_observables(m)

    # Initialize pseudo ZZ and DD matrices
    ZZ_pseudo = zeros(_n_pseudo, _n_states)
    DD_pseudo = zeros(_n_pseudo)

    no_integ_inds = inds_states_no_integ_series(m)
    if get_setting(m, :add_laborproductivity_measurement)
        # Construct pseudo-obs from integrated states first
        ZZ_pseudo[pseudo[:laborproductivity], endo[:y_t]] = 1.
        ZZ_pseudo[pseudo[:laborproductivity], endo[:L_t]] = -1.
        # ZZ_pseudo[pseudo[:laborproductivity], endo_addl[:cum_z_t]] = 1.
        DD_pseudo[pseudo[:laborproductivity]] = 100. * log(m[:ystar] / m[:Lstar])

        # Remove integrated states (e.g. states w/unit roots)
        # RRR and CCC aren't used, so we don't do anything with them
        TTT = @view TTT[no_integ_inds, no_integ_inds]
    end

    # Compute TTT^10, used for Expected10YearRateGap, Expected10YearRate, and Expected10YearNaturalRate
    TTT10 = (1/40)*((UniformScaling(1.) - TTT)\(TTT - TTT^41))

    ##########################################################
    ## PSEUDO-OBSERVABLE EQUATIONS
    ##########################################################

    ## Output
    ZZ_pseudo[pseudo[:y_t],endo[:y_t]] = 1.

    ## Flexible Output
    ZZ_pseudo[pseudo[:y_f_t],endo[:y_f_t]] = 1.

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

    ## Hours
    ZZ_pseudo[pseudo[:Hours],endo[:L_t]] = 1.

    ## Flexible Hours
    ZZ_pseudo[pseudo[:FlexibleHours],endo[:L_f_t]] = 1.

    ## z_t
    ZZ_pseudo[pseudo[:z_t], endo[:z_t]] = 1.

    ## Expected 10-Year Rate Gap
    ZZ_pseudo[pseudo[:Expected10YearRateGap], no_integ_inds] = TTT10[endo[:R_t], :] - TTT10[endo[:r_f_t], :] - TTT10[endo[:Eπ_t], :]

    ## Nominal FFR
    ZZ_pseudo[pseudo[:NominalFFR], endo[:R_t]] = 1.
    DD_pseudo[pseudo[:NominalFFR]] = m[:Rstarn]

    ## Expected 10-Year Interest Rate
    ZZ_pseudo[pseudo[:Expected10YearRate], no_integ_inds] = TTT10[endo[:R_t], :]
    DD_pseudo[pseudo[:Expected10YearRate]]    = m[:Rstarn]

    ## Expected 10-Year Natural Rate
    ZZ_pseudo[pseudo[:Expected10YearNaturalRate], no_integ_inds] = TTT10[endo[:r_f_t], :] + TTT10[endo[:Eπ_t], :]
    DD_pseudo[pseudo[:Expected10YearNaturalRate]]    = m[:Rstarn]

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
    ZZ_pseudo[pseudo[:NominalWageGrowth],endo[:w_t]] = 1.
    ZZ_pseudo[pseudo[:NominalWageGrowth],endo_addl[:w_t1]] = -1.
    ZZ_pseudo[pseudo[:NominalWageGrowth],endo[:z_t]] = 1.
    ZZ_pseudo[pseudo[:NominalWageGrowth],endo[:π_t]] = 1.
    DD_pseudo[pseudo[:NominalWageGrowth]]            = 100*(m[:π_star]-1) + 100*(exp(m[:z_star])-1)

    ## labor share
    if haskey(m.settings, :add_laborshare_measurement)
        if get_setting(m, :add_laborshare_measurement)
            ZZ_pseudo[pseudo[:laborshare_t], endo[:w_t]] = 1.
            ZZ_pseudo[pseudo[:laborshare_t], endo[:L_t]] = 1.
            ZZ_pseudo[pseudo[:laborshare_t], endo[:y_t]] = -1.
            DD_pseudo[pseudo[:laborshare_t]] = 100. * log(m[:wstar] * m[:Lstar] / m[:ystar])
        end
    end


    ## Labor Productivity Growth, no measurement error
    if haskey(m.settings, :add_laborproductivitygrowth_nome_measurement)
        if get_setting(m, :add_laborproductivitygrowth_nome_measurement)
            ZZ_pseudo[pseudo[:LaborProductivityGrowthNoME], endo[:y_t]]       = 1.
            ZZ_pseudo[pseudo[:LaborProductivityGrowthNoME], endo_addl[:y_t1]] = -1.
            ZZ_pseudo[pseudo[:LaborProductivityGrowthNoME], endo[:z_t]]       = 1.
            ZZ_pseudo[pseudo[:LaborProductivityGrowthNoME], endo[:L_t]]       = -1
            ZZ_pseudo[pseudo[:LaborProductivityGrowthNoME], endo_addl[:L_t1]] = 1.
            DD_pseudo[pseudo[:LaborProductivityGrowthNoME]]                   = 100*(exp(m[:z_star]) - 1)
        end
    end

    ## Fundamental inflation related pseudo-obs
    if subspec(m) in ["ss13", "ss14", "ss15", "ss16", "ss17", "ss18", "ss19"]
        # Compute coefficient on Sinf
        betabar = exp((1-m[:σ_c] ) * m[:z_star]) * m[:β]
        if subspec in ["ss21", "ss22", "ss25", "ss26"] && reg == 2
            κ = ((1 - m[:ζ_p_r2]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))*
                 (1 - m[:ζ_p_r2]))/(m[:ζ_p_r2]*((m[:Φ]- 1)*m[:ϵ_p] + 1))/
                 (1 + m[:ι_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        else
            κ = ((1 - m[:ζ_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))*
                 (1 - m[:ζ_p]))/(m[:ζ_p]*((m[:Φ]- 1)*m[:ϵ_p] + 1))/
                 (1 + m[:ι_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        end
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

    return PseudoMeasurement(ZZ_pseudo, DD_pseudo)
end


function pseudo_measurement(m::Model1002{T},
                            TTTs::Vector{Matrix{T}},
                            RRRs::Vector{Matrix{T}},
                            CCCs::Vector{Vector{T}}) where {T<:AbstractFloat}

    endo      = m.endogenous_states
    endo_addl = m.endogenous_states_augmented
    pseudo    = m.pseudo_observables

    _n_states = n_states_augmented(m)
    _n_pseudo = n_pseudo_observables(m)

    n_reg = length(TTTs)

    ZZ_pseudos = Vector{Matrix{Float64}}(undef, n_reg)
    DD_pseudos = Vector{Vector{Float64}}(undef, n_reg)

    for reg in 1:n_reg
        # Initialize pseudo ZZ and DD matrices
        ZZ_pseudos[reg] = zeros(_n_pseudo, _n_states)
        DD_pseudos[reg] = zeros(_n_pseudo)

        no_integ_inds = inds_states_no_integ_series(m)
        if get_setting(m, :add_laborproductivity_measurement)
            # Construct pseudo-obs from integrated states first
            ZZ_pseudos[reg][pseudo[:laborproductivity], endo[:y_t]] = 1.
            ZZ_pseudos[reg][pseudo[:laborproductivity], endo[:L_t]] = -1.
            ZZ_pseudos[reg][pseudo[:laborproductivity], endo_addl[:cum_z_t]] = 1.
            DD_pseudos[reg][pseudo[:laborproductivity]] = 100. * log(m[:ystar] / m[:Lstar])

            # Remove integrated states (e.g. states w/unit roots)
            # RRR and CCC aren't used, so we don't do anything with them
            TTT = @view TTTs[reg][no_integ_inds, no_integ_inds]
        end

        # Compute TTT^10, used for Expected10YearRateGap, Expected10YearRate, and Expected10YearNaturalRate
        TTT10 = (1/40)*((UniformScaling(1.) - TTTs[reg])\(TTTs[reg] - TTTs[reg]^41))

        ##########################################################
        ## PSEUDO-OBSERVABLE EQUATIONS
        ##########################################################

        ## Output
        ZZ_pseudos[reg][pseudo[:y_t],endo[:y_t]] = 1.

        ## Flexible Output
        ZZ_pseudos[reg][pseudo[:y_f_t],endo[:y_f_t]] = 1.

        ## Natural Rate
        ZZ_pseudos[reg][pseudo[:NaturalRate],endo[:r_f_t]] = 1.
        DD_pseudos[reg][pseudo[:NaturalRate]]              = 100.0*(m[:rstar]-1.0)

        ## π_t
        ZZ_pseudos[reg][pseudo[:π_t],endo[:π_t]] = 1.
        DD_pseudos[reg][pseudo[:π_t]]            = 100*(m[:π_star]-1);

        ## Output Gap
        ZZ_pseudos[reg][pseudo[:OutputGap],endo[:y_t]] = 1;
        ZZ_pseudos[reg][pseudo[:OutputGap],endo[:y_f_t]] = -1;

        ## Ex Ante Real Rate
        ZZ_pseudos[reg][pseudo[:ExAnteRealRate],endo[:R_t]]  = 1;
        ZZ_pseudos[reg][pseudo[:ExAnteRealRate],endo[:Eπ_t]] = -1;
        DD_pseudos[reg][pseudo[:ExAnteRealRate]]             = m[:Rstarn] - 100*(m[:π_star]-1);

        ## Long Run Inflation
        ZZ_pseudos[reg][pseudo[:LongRunInflation],endo[:π_star_t]] = 1.
        DD_pseudos[reg][pseudo[:LongRunInflation]]                 = 100. *(m[:π_star]-1.)

        ## Marginal Cost
        ZZ_pseudos[reg][pseudo[:MarginalCost],endo[:mc_t]] = 1.

        ## Wages
        ZZ_pseudos[reg][pseudo[:Wages],endo[:w_t]] = 1.

        ## Flexible Wages
        ZZ_pseudos[reg][pseudo[:FlexibleWages],endo[:w_f_t]] = 1.

        ## Hours
        ZZ_pseudos[reg][pseudo[:Hours],endo[:L_t]] = 1.

        ## Flexible Hours
        ZZ_pseudos[reg][pseudo[:FlexibleHours],endo[:L_f_t]] = 1.

        ## z_t
        ZZ_pseudos[reg][pseudo[:z_t], endo[:z_t]] = 1.

        ## Expected 10-Year Rate Gap
        ZZ_pseudos[reg][pseudo[:Expected10YearRateGap], no_integ_inds] = TTT10[endo[:R_t], :] - TTT10[endo[:r_f_t], :] - TTT10[endo[:Eπ_t], :]

        ## Nominal FFR
        ZZ_pseudos[reg][pseudo[:NominalFFR], endo[:R_t]] = 1.
        DD_pseudos[reg][pseudo[:NominalFFR]] = m[:Rstarn]

        ## Expected 10-Year Interest Rate
        ZZ_pseudos[reg][pseudo[:Expected10YearRate], no_integ_inds] = TTT10[endo[:R_t], :]
        DD_pseudos[reg][pseudo[:Expected10YearRate]]    = m[:Rstarn]

        ## Expected 10-Year Natural Rate
        ZZ_pseudos[reg][pseudo[:Expected10YearNaturalRate], no_integ_inds] = TTT10[endo[:r_f_t], :] + TTT10[endo[:Eπ_t], :]
        DD_pseudos[reg][pseudo[:Expected10YearNaturalRate]]    = m[:Rstarn]

        ## Expected Nominal Natural Rate
        ZZ_pseudos[reg][pseudo[:ExpectedNominalNaturalRate], endo[:r_f_t]] = 1.
        ZZ_pseudos[reg][pseudo[:ExpectedNominalNaturalRate], endo[:Eπ_t]]  = 1.
        DD_pseudos[reg][pseudo[:ExpectedNominalNaturalRate]]               = m[:Rstarn]

        ## Nominal Rate Gap
        ZZ_pseudos[reg][pseudo[:NominalRateGap], endo[:R_t]]   = 1.
        ZZ_pseudos[reg][pseudo[:NominalRateGap], endo[:r_f_t]] = -1.
        ZZ_pseudos[reg][pseudo[:NominalRateGap], endo[:Eπ_t]]  = -1.

        ## Labor Productivity Growth
        ZZ_pseudos[reg][pseudo[:LaborProductivityGrowth], endo[:y_t]]           = 1.
        ZZ_pseudos[reg][pseudo[:LaborProductivityGrowth], endo_addl[:y_t1]]     = -1.
        ZZ_pseudos[reg][pseudo[:LaborProductivityGrowth], endo[:z_t]]           = 1.
        ZZ_pseudos[reg][pseudo[:LaborProductivityGrowth], endo_addl[:e_gdp_t]]  = 1.
        ZZ_pseudos[reg][pseudo[:LaborProductivityGrowth], endo_addl[:e_gdp_t1]] = -m[:me_level]
        ZZ_pseudos[reg][pseudo[:LaborProductivityGrowth], endo[:L_t]]           = -1
        ZZ_pseudos[reg][pseudo[:LaborProductivityGrowth], endo_addl[:L_t1]]     = 1.
        DD_pseudos[reg][pseudo[:LaborProductivityGrowth]]                       = 100*(exp(m[:z_star]) - 1)

        ## u_t
        ZZ_pseudos[reg][pseudo[:u_t], endo[:u_t]] = 1.

        ## Nominal Wage Growth
        ZZ_pseudos[reg][pseudo[:NominalWageGrowth],endo[:w_t]] = 1.
        ZZ_pseudos[reg][pseudo[:NominalWageGrowth],endo_addl[:w_t1]] = -1.
        ZZ_pseudos[reg][pseudo[:NominalWageGrowth],endo[:π_t]] = 1.
        DD_pseudos[reg][pseudo[:NominalWageGrowth]]            = 100*(m[:π_star]-1);

        ## labor share
        if haskey(m.settings, :add_laborshare_measurement)
            if get_setting(m, :add_laborshare_measurement)
                ZZ_pseudos[reg][pseudo[:laborshare_t], endo[:w_t]] = 1.
                ZZ_pseudos[reg][pseudo[:laborshare_t], endo[:L_t]] = 1.
                ZZ_pseudos[reg][pseudo[:laborshare_t], endo[:y_t]] = -1.
                DD_pseudos[reg][pseudo[:laborshare_t]] = 100. * log(m[:wstar] * m[:Lstar] / m[:ystar])
            end
        end

<<<<<<< Updated upstream
        ## Labor Productivity Growth, no measurement error
        if haskey(m.settings, :add_laborproductivitygrowth_nome_measurement)
            if get_setting(m, :add_laborproductivitygrowth_nome_measurement)
                ZZ_pseudos[reg][pseudo[:LaborProductivityGrowthNoME], endo[:y_t]]       = 1.
                ZZ_pseudos[reg][pseudo[:LaborProductivityGrowthNoME], endo_addl[:y_t1]] = -1.
                ZZ_pseudos[reg][pseudo[:LaborProductivityGrowthNoME], endo[:z_t]]       = 1.
                ZZ_pseudos[reg][pseudo[:LaborProductivityGrowthNoME], endo[:L_t]]       = -1
                ZZ_pseudos[reg][pseudo[:LaborProductivityGrowthNoME], endo_addl[:L_t1]] = 1.
                DD_pseudos[reg][pseudo[:LaborProductivityGrowthNoME]]                   = 100*(exp(m[:z_star]) - 1)
            end
        end

        ## Fundameantal inflation related pseudo-obs
        if subspec(m) in ["ss13", "ss14", "ss15", "ss16", "ss17", "ss18", "ss19", "ss20"] #,
                          # "ss21", "ss22", "ss23","ss24", "ss25", "ss26"]
            # Compute coefficient on Sinf
            betabar = exp((1-m[:σ_c] ) * m[:z_star]) * m[:β]
            if subspec(m) in ["ss21", "ss22", "ss25", "ss26", "ss28", "ss29", "ss41", "ss42"] && reg == 2
                κ = ((1 - m[:ζ_p_r2]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))*
                     (1 - m[:ζ_p_r2]))/(m[:ζ_p_r2]*((m[:Φ]- 1)*m[:ϵ_p] + 1))/
                (1 + m[:ι_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
            else
                κ = ((1 - m[:ζ_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))*
                     (1 - m[:ζ_p]))/(m[:ζ_p]*((m[:Φ]- 1)*m[:ϵ_p] + 1))/
                (1 + m[:ι_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
            end

            κcoef = κ * (1 + m[:ι_p] * betabar)

            ZZ_pseudos[reg][pseudo[:Sinf_t], endo_addl[:Sinf_t]] = 1.
            ZZ_pseudos[reg][pseudo[:Sinf_w_coef_t], endo_addl[:Sinf_t]] = κcoef
            DD_pseudos[reg][pseudo[:ι_p]] = m[:ι_p]
            ZZ_pseudos[reg][pseudo[:πtil_t], endo_addl[:πtil_t]] = 1.
            DD_pseudos[reg][pseudo[:πtil_t]] = 100 * (m[:π_star] - 1)
            ZZ_pseudos[reg][pseudo[:e_tfp_t], endo_addl[:e_tfp_t]] = 1.
            if subspec(m) in ["ss14", "ss15", "ss16", "ss18", "ss19"]
                ZZ_pseudos[reg][pseudo[:e_tfp_t1], endo_addl[:e_tfp_t1]] = 1.
            end
        end

        ## Exogenous processes
        if subspec(m) == "ss12"
            to_add = [:g_t, :b_t, :μ_t, :z_t, :λ_f_t, :λ_w_t, :rm_t, :σ_ω_t, :μ_e_t,
                      :γ_t, :π_star_t]
            to_add_addl = [:e_lr_t, :e_tfp_t, :e_gdpdef_t, :e_corepce_t, :e_gdp_t, :e_gdi_t]
            for i in to_add
                ZZ_pseudos[reg][pseudo[i], endo[i]] = 1.
            end
            for i in to_add_addl
                ZZ_pseudos[reg][pseudo[i], endo_addl[i]] = 1.
            end
        end
    end
    return [PseudoMeasurement(ZZ_pseudos[1], DD_pseudos[1]), PseudoMeasurement(ZZ_pseudos[2], DD_pseudos[2])]
end
