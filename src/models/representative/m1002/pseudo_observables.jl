function init_pseudo_observable_mappings!(m::Model1002)

    pseudo_names = [:y_t, :y_f_t, :NaturalRate, :π_t, :OutputGap, :ExAnteRealRate, :LongRunInflation,
                    :MarginalCost, :Wages, :FlexibleWages, :Hours, :FlexibleHours, :z_t,
                    :Expected10YearRateGap, :NominalFFR, :Expected10YearRate,
                    :Expected10YearNaturalRate,
                    :ExpectedNominalNaturalRate, :NominalRateGap, :LaborProductivityGrowth, :u_t]
                    # :i_f_t, :R_t, :c_f_t, :qk_f_t, :k_f_t, :kbar_f_t, :u_f_t, :rk_f_t, :w_f_t, :L_f_t, :rktil_f_t, :n_f_t, :c_t,
                    # :b_t, :r_f_t]

    if subspec(m) == "ss12"
        to_add = [:g_t, :b_t, :μ_t, :λ_f_t, :λ_w_t, :rm_t, :σ_ω_t, :μ_e_t,
                  :γ_t, :π_star_t, :e_lr_t, :e_tfp_t, :e_gdpdef_t, :e_corepce_t, :e_gdp_t, :e_gdi_t]
        pseudo_names = vcat(pseudo_names, to_add)
    end

    if subspec(m) in ["ss13", "ss14", "ss15", "ss16", "ss17", "ss18", "ss19"]
        push!(pseudo_names, :Sinf_t, :Sinf_w_coef_t, :ι_p, :πtil_t, :πtil_t1, :e_tfp_t)
        if subspec(m) in ["ss14", "ss15", "ss16", "ss18", "ss19"]
            push!(pseudo_names, :e_tfp_t1)
        end
    end

    if haskey(get_settings(m), :add_covid_pseudoobs)
        if get_setting(m, :add_covid_pseudoobs) && subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63", "ss64", "ss65", "ss66", "ss67", "ss68", "ss69", "ss70", "ss71", "ss72", "ss73", "ss74", "ss75", "ss76", "ss77", "ss78", "ss79", "ss80", "ss81", "ss82", "ss83", "ss84", "ss85"]
            push!(pseudo_names, :ziid, :varphiiid, :biidc)
        end
    end

    if haskey(get_settings(m), :add_pseudo_gdp)
        if get_setting(m, :add_pseudo_gdp) && subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63", "ss64", "ss65", "ss66", "ss67", "ss68", "ss69", "ss70", "ss71", "ss72", "ss73", "ss74", "ss75", "ss76", "ss77", "ss78", "ss79", "ss80", "ss81", "ss82", "ss83", "ss84", "ss85"]
            push!(pseudo_names, :PseudoGDP)
        end
    end

    if haskey(get_settings(m), :add_pseudo_corepce)
        if get_setting(m, :add_pseudo_corepce) && subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63", "ss64", "ss65", "ss66", "ss67", "ss68", "ss69", "ss70", "ss71", "ss72", "ss73", "ss74", "ss75", "ss76", "ss77", "ss78", "ss79", "ss80", "ss81", "ss82", "ss83", "ss84", "ss85"]
            push!(pseudo_names, :PseudoCorePCE)
        end
    end

    if haskey(get_settings(m), :add_NominalWageGrowth)
        if get_setting(m, :add_NominalWageGrowth)
            push!(pseudo_names, :NominalWageGrowth)
        end
    end
    if haskey(get_settings(m), :add_ztil)
        if get_setting(m, :add_ztil)
            push!(pseudo_names, :ztil)
        end
    end
    if haskey(get_settings(m), :add_zp)
        if get_setting(m, :add_zp)
            push!(pseudo_names, :zp)
        end
    end
    if haskey(get_settings(m), :add_pgap)
        if get_setting(m, :add_pgap)
            push!(pseudo_names, :pgap)
        end
    end
    if haskey(get_settings(m), :add_ygap)
        if get_setting(m, :add_ygap)
            push!(pseudo_names, :ygap)
        end
    end
    if haskey(get_settings(m), :add_altpolicy_pgap)
        if get_setting(m, :add_altpolicy_pgap)
            push!(pseudo_names, :pgap)
        end
    end
    if haskey(get_settings(m), :add_altpolicy_ygap)
        if get_setting(m, :add_altpolicy_ygap)
            push!(pseudo_names, :ygap)
        end
    end

    if haskey(get_settings(m), :add_urhat)
        if get_setting(m, :add_urhat)
            push!(pseudo_names, :urhat)
        end
    end

    if haskey(get_settings(m), :add_rw)
        if get_setting(m, :add_rw)
            push!(pseudo_names, :rw)
            push!(pseudo_names, :Rref)
        end
    end

    if haskey(get_settings(m), :add_laborshare_measurement)
        if get_setting(m, :add_laborshare_measurement)
            push!(pseudo_names, :laborshare_t)
        end
    end

    if haskey(get_settings(m), :add_laborproductivity_measurement)
        if get_setting(m, :add_laborproductivity_measurement)
            push!(pseudo_names, :laborproductivity)
        end
    end

    if haskey(get_settings(m), :add_nominalgdp_level)
        if get_setting(m, :add_nominalgdp_level)
            push!(pseudo_names, :NominalGDPLevel)
        end
    end

    if haskey(get_settings(m), :add_nominalgdp_growth)
        if get_setting(m, :add_nominalgdp_growth)
            push!(pseudo_names, :NominalGDPGrowth)
        end
    end

    if haskey(get_settings(m), :add_flexible_price_growth)
        if get_setting(m, :add_flexible_price_growth)
            push!(pseudo_names, :FlexibleGDPGrowth, :FlexibleConsumptionGrowth, :FlexibleInvestmentGrowth)
        end
    end

    if haskey(get_settings(m), :add_cumulative)
        if get_setting(m, :add_cumulative)
            push!(pseudo_names, :AccumOutputGap, :GDPLevel,
                  :TechnologyLevel, :FlexibleGDPLevel, :ConsumptionLevel,
                  :FlexibleConsumptionLevel, :InvestmentLevel,
                  :FlexibleInvestmentLevel)
        end
    end

    if haskey(get_settings(m), :add_laborproductivitygrowth_nome_measurement)
        if get_setting(m, :add_laborproductivitygrowth_nome_measurement)
            push!(pseudo_names, :LaborProductivityGrowthNoME)
        end
    end

    if haskey(get_settings(m), :add_Epi_t_measurement)
        if get_setting(m, :add_Epi_t_measurement)
            push!(pseudo_names, :Epi_t)
        end
    end

    # Create PseudoObservable objects
    pseudo = OrderedDict{Symbol,PseudoObservable}()
    for k in pseudo_names
        pseudo[k] = PseudoObservable(k)
    end

    # Fill in names and reverse transforms
    pseudo[:y_t].name = "Output Growth"
    pseudo[:y_t].longname = "Output Growth Per Capita"

    pseudo[:y_f_t].name = "Flexible Output Growth"
    pseudo[:y_f_t].longname = "Output that would prevail in a flexible-price economy."

    pseudo[:NaturalRate].name = "Real Natural Rate"
    pseudo[:NaturalRate].longname = "The real interest rate that would prevail in a flexible-price economy."
    pseudo[:NaturalRate].rev_transform = quartertoannual

    pseudo[:π_t].name = "Inflation"
    pseudo[:π_t].longname = "Inflation"
    pseudo[:π_t].rev_transform = quartertoannual

    pseudo[:OutputGap].name = "Output Gap"
    pseudo[:OutputGap].longname = "Output Gap"

    pseudo[:ExAnteRealRate].name = "Ex Ante Real Rate"
    pseudo[:ExAnteRealRate].longname = "Ex Ante Real Rate"
    pseudo[:ExAnteRealRate].rev_transform = quartertoannual

    pseudo[:LongRunInflation].name = "Long Run Inflation"
    pseudo[:LongRunInflation].longname = "Long Run Inflation"
    pseudo[:LongRunInflation].rev_transform = quartertoannual

    pseudo[:MarginalCost].name = "Marginal Cost"
    pseudo[:MarginalCost].longname = "Marginal Cost"

    pseudo[:Wages].name = "Wages"
    pseudo[:Wages].longname = "Wages"

    pseudo[:FlexibleWages].name = "Flexible Wages"
    pseudo[:FlexibleWages].longname = "Wages that would prevail in a flexible-wage economy"

    pseudo[:Hours].name = "Hours"
    pseudo[:Hours].longname = "Hours"

    pseudo[:FlexibleHours].name     = "Flexible Hours"
    pseudo[:FlexibleHours].longname = "Flexible Hours"

    pseudo[:z_t].name     = "z_t (Technology Growth minus Steady State Growth)"
    pseudo[:z_t].longname = "z_t (Technology Growth minus Steady State Growth)"

    pseudo[:Expected10YearRateGap].name     = "Expected 10-Year Rate Gap"
    pseudo[:Expected10YearRateGap].longname = "Expected 10-Year Rate Gap"
    pseudo[:Expected10YearRateGap].rev_transform = quartertoannual

    pseudo[:NominalFFR].name     = "Nominal FFR"
    pseudo[:NominalFFR].longname = "Nominal FFR at an annual rate"
    pseudo[:NominalFFR].rev_transform = quartertoannual

    pseudo[:Expected10YearRate].name     = "Expected 10-Year Rate"
    pseudo[:Expected10YearRate].longname = "Expected 10-Year Interest Rate"
    pseudo[:Expected10YearRate].rev_transform = quartertoannual

    pseudo[:Expected10YearNaturalRate].name     = "Expected 10-Year Natural Rate"
    pseudo[:Expected10YearNaturalRate].longname = "Expected 10-Year Natural Rate of Interest"
    pseudo[:Expected10YearNaturalRate].rev_transform = quartertoannual

    pseudo[:ExpectedNominalNaturalRate].name     = "Expected Nominal Natural Rate"
    pseudo[:ExpectedNominalNaturalRate].longname = "Natural Rate + Expected Inflation"
    pseudo[:ExpectedNominalNaturalRate].rev_transform = quartertoannual

    pseudo[:NominalRateGap].name     = "Nominal Rate Gap"
    pseudo[:NominalRateGap].longname = "Nominal FFR - Nominal Natural Rate"
    pseudo[:NominalRateGap].rev_transform = quartertoannual

    pseudo[:LaborProductivityGrowth].name     = "Labor Productivity Growth"
    pseudo[:LaborProductivityGrowth].longname = "Labor Productivity Growth Rate"
    pseudo[:LaborProductivityGrowth].rev_transform = quartertoannual

    pseudo[:u_t].name     = "u_t"
    pseudo[:u_t].longname = "u_t"

    # pseudo[:i_f_t].name     = "i_f_t"
    # pseudo[:i_f_t].longname = "i_f_t"

    # pseudo[:R_t].name     = "R_t"
    # pseudo[:R_t].longname = "R_t"

    # pseudo[:c_f_t].name     = "c_f_t"
    # pseudo[:c_f_t].longname = "c_f_t"

    # pseudo[:r_f_t].name     = "r_f_t"
    # pseudo[:r_f_t].longname = "r_f_t"

    # pseudo[:qk_f_t].name     = "qk_f_t"
    # pseudo[:qk_f_t].longname = "qk_f_t"

    # pseudo[:k_f_t].name     = "k_f_t"
    # pseudo[:k_f_t].longname = "k_f_t"

    # pseudo[:kbar_f_t].name     = "kbar_f_t"
    # pseudo[:kbar_f_t].longname = "kbar_f_t"

    # pseudo[:u_f_t].name     = "u_f_t"
    # pseudo[:u_f_t].longname = "u_f_t"

    # pseudo[:rk_f_t].name     = "rk_f_t"
    # pseudo[:rk_f_t].longname = "rk_f_t"

    # pseudo[:w_f_t].name     = "w_f_t"
    # pseudo[:w_f_t].longname = "w_f_t"

    # pseudo[:L_f_t].name     = "L_f_t"
    # pseudo[:L_f_t].longname = "L_f_t"


    # pseudo[:rktil_f_t].name     = "rktil_f_t"
    # pseudo[:rktil_f_t].longname = "rktil_f_t"

    # pseudo[:b_t].name     = "b_t"
    # pseudo[:b_t].longname = "b_t"

    # pseudo[:n_f_t].name     = "n_f_t"
    # pseudo[:n_f_t].longname = "n_f_t"

    if subspec(m) in ["ss13", "ss14", "ss15", "ss16", "ss17", "ss18", "ss19"]
        pseudo[:Sinf_t].name     = "Sinf_t"
        pseudo[:Sinf_t].longname = "Sinf_t, PDV of Emc_t"
        pseudo[:Sinf_w_coef_t].name     = "Sinf_w_coef_t"
        pseudo[:Sinf_w_coef_t].longname = "Sinf_w_coef_t, PDV of Emc_t multiplied by coefficient"
        pseudo[:ι_p].name = "iota_p"
        pseudo[:ι_p].longname = "iota_p"
        pseudo[:πtil_t].name     = "pitil_t"
        pseudo[:πtil_t].longname = "Fundamental Inflation"
        pseudo[:πtil_t1].name     = "pitil_t1"
        pseudo[:πtil_t1].longname = "Fundamental Inflation Lag 1"
        pseudo[:e_tfp_t].name = "e_tfp_t"
        pseudo[:e_tfp_t].longname = "e_tfp_t"
        if subspec(m) in ["ss14", "ss15", "ss16", "ss18", "ss19"]
            pseudo[:e_tfp_t1].name = "e_tfp_t1"
            pseudo[:e_tfp_t1].longname = "e_tfp_t1"
        end
    end

    # Other exogenous processes
    if subspec(m) == "ss12"
        for i in to_add
            pseudo[i].name = string(DSGE.detexify(i))
            pseudo[i].longname = string(DSGE.detexify(i))
        end
    end

    if haskey(get_settings(m), :add_NominalWageGrowth)
        if get_setting(m, :add_NominalWageGrowth)
            pseudo[:NominalWageGrowth].name = "Nominal Wage Growth"
            pseudo[:NominalWageGrowth].longname = "Nominal Wage Growth"
        end
    end

    if haskey(get_settings(m), :add_laborshare_measurement)
        if get_setting(m, :add_laborshare_measurement)
            pseudo[:laborshare_t].name     = "Log Labor Share"
            pseudo[:laborshare_t].longname = "Log Labor Share"
        end
    end

    if haskey(get_settings(m), :add_nominalgdp_level)
        if get_setting(m, :add_nominalgdp_level)
            pseudo[:NominalGDPLevel].name     = "Nominal GDP Level"
            pseudo[:NominalGDPLevel].longname = "Nominal GDP Level"
        end
    end

    if haskey(get_settings(m), :add_nominalgdp_growth)
        if get_setting(m, :add_nominalgdp_growth)
            pseudo[:NominalGDPGrowth].name     = "Nominal GDP Growth"
            pseudo[:NominalGDPGrowth].longname = "Nominal GDP Growth"
        end
    end

    if haskey(get_settings(m), :add_cumulative)
        if get_setting(m, :add_cumulative)
            pseudo[:AccumOutputGap].name     = "Accumulated Output Gap"
            pseudo[:AccumOutputGap].longname = "Accumulated Output Gap"
            pseudo[:GDPLevel].name     = "GDP Level"
            pseudo[:GDPLevel].longname = "GDP Level"
            pseudo[:TechnologyLevel].name     = "Technology (z) Level"
            pseudo[:TechnologyLevel].longname = "Technology (z) Level"
            pseudo[:FlexibleGDPLevel].name     = "Flexible GDP Level"
            pseudo[:FlexibleGDPLevel].longname = "Flexible GDP Level"
            pseudo[:ConsumptionLevel].name     = "Consumption Level"
            pseudo[:ConsumptionLevel].longname = "Consumption Level"
            pseudo[:FlexibleConsumptionLevel].name     = "Flexible Consumption Level"
            pseudo[:FlexibleConsumptionLevel].longname = "Flexible Consumption Level"
            pseudo[:InvestmentLevel].name     = "Investment Level"
            pseudo[:InvestmentLevel].longname = "Investment Level"
            pseudo[:FlexibleInvestmentLevel].name     = "Flexible Investment Level"
            pseudo[:FlexibleInvestmentLevel].longname = "Flexible Investment Level"
        end
    end

    if haskey(get_settings(m), :add_flexible_price_growth)
        if get_setting(m, :add_flexible_price_growth)
            pseudo[:FlexibleGDPGrowth].name     = "Flexible GDP Level"
            pseudo[:FlexibleGDPGrowth].longname = "Flexible GDP Growth"
            pseudo[:FlexibleConsumptionGrowth].name     = "Flexible Consumption Growth"
            pseudo[:FlexibleConsumptionGrowth].longname = "Flexible Consumption Growth"
            pseudo[:FlexibleInvestmentGrowth].name     = "Flexible Investment Growth"
            pseudo[:FlexibleInvestmentGrowth].longname = "Flexible Investment Growth"
        end
    end

    if haskey(get_settings(m), :add_ztil)
        if get_setting(m, :add_ztil)
            pseudo[:ztil].name     = "ztil"
            pseudo[:ztil].longname = "ztil"
        end
    end

    if haskey(get_settings(m), :add_zp)
        if get_setting(m, :add_zp)
            pseudo[:zp].name     = "zp"
            pseudo[:zp].longname = "zp"
        end
    end

    if haskey(get_settings(m), :add_altpolicy_pgap)
        if get_setting(m, :add_altpolicy_pgap)
            pseudo[:pgap].name     = "pgap"
            pseudo[:pgap].longname = "pgap"
        end
    end

    if haskey(get_settings(m), :add_altpolicy_ygap)
        if get_setting(m, :add_altpolicy_ygap)
            pseudo[:ygap].name     = "ygap"
            pseudo[:ygap].longname = "ygap"
        end
    end

    if haskey(get_settings(m), :add_urhat)
        if get_setting(m, :add_urhat)
            pseudo[:urhat].name     = "urhat"
            pseudo[:urhat].longname = "urhat"
        end
    end

    if haskey(get_settings(m), :add_rw)
        if get_setting(m, :add_rw)
            pseudo[:rw].name     = "rw"
            pseudo[:rw].longname = "rw"
            pseudo[:Rref].name     = "Rref"
            pseudo[:Rref].longname = "Rref"
        end
    end

    if haskey(get_settings(m), :add_covid_pseudoobs)
        if get_setting(m, :add_covid_pseudoobs) && subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63", "ss64", "ss65", "ss66", "ss67", "ss68", "ss69", "ss70", "ss71", "ss72", "ss73", "ss74", "ss75", "ss76", "ss77", "ss78", "ss79", "ss80", "ss81", "ss82", "ss83", "ss84", "ss85"]
            pseudo[:ziid].name     = "ziid"
            pseudo[:ziid].longname = "ziid"
            pseudo[:biidc].name     = "biidc"
            pseudo[:biidc].longname = "biidc"
            pseudo[:varphiiid].name     = "varphiiid"
            pseudo[:varphiiid].longname = "varphiiid"
        end
    end

    if haskey(get_settings(m), :add_pseudo_gdp)
        if get_setting(m, :add_pseudo_gdp) && subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63", "ss64", "ss65", "ss66", "ss67", "ss68", "ss69", "ss70", "ss71", "ss72", "ss73", "ss74", "ss75", "ss76", "ss77", "ss78", "ss79", "ss80", "ss81", "ss82", "ss83", "ss84", "ss85"]
            pseudo[:PseudoGDP].name     = "GDP Growth Pseudo-observable"
            pseudo[:PseudoGDP].longname = "GDP Growth Pseudo-observable"
            pseudo[:PseudoGDP].rev_transform = loggrowthtopct_annualized_percapita
        end
    end

    if haskey(get_settings(m), :add_pseudo_corepce)
        if get_setting(m, :add_pseudo_corepce) && subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63", "ss64", "ss65", "ss66", "ss67", "ss68", "ss69", "ss70", "ss71", "ss72", "ss73", "ss74", "ss75", "ss76", "ss77", "ss78", "ss79", "ss80", "ss81", "ss82", "ss83", "ss84", "ss85"]
            pseudo[:PseudoCorePCE].name     = "Core PCE Pseudo-observable"
            pseudo[:PseudoCorePCE].longname = "Core PCE Pseudo-observable"
            pseudo[:PseudoCorePCE].rev_transform = loggrowthtopct_annualized
        end
    end

    if haskey(get_settings(m), :add_laborproductivity_measurement)
        if get_setting(m, :add_laborproductivity_measurement)
            pseudo[:laborproductivity].name     = "Log Labor Productivity"
            pseudo[:laborproductivity].longname = "Log Labor Productivity"
        end
    end

    if haskey(get_settings(m), :add_laborproductivitygrowth_nome_measurement)
        if get_setting(m, :add_laborproductivitygrowth_nome_measurement)
            pseudo[:LaborProductivityGrowthNoME].name     = "Labor Productivity Growth (No ME)"
            pseudo[:LaborProductivityGrowthNoME].longname = "Labor Productivity Growth (No ME)"
        end
    end

    if haskey(get_settings(m), :add_Epi_t_measurement)
        if get_setting(m, :add_Epi_t_measurement)
            pseudo[:Epi_t].name     = "Short-term Inflation Expectations"
            pseudo[:Epi_t].longname = "Short-term Inflation Expectations"
        end
    end

    # Needed to implement pseudo-measurement equation correctly
    m <= Setting(:forward_looking_pseudo_observables, [:Expected10YearRateGap, :Expected10YearRate, :Expected10YearNaturalRate])

    # Add to model object
    m.pseudo_observable_mappings = pseudo
end
