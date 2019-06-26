function init_observable_mappings!(m::RealBondMkup)

    observables = OrderedDict{Symbol,Observable}()
    population_mnemonic = get(get_setting(m, :population_mnemonic))

    ############################################################################
    ## 1. GDP
    ############################################################################
    gdp_fwd_transform = function (levels)
        # FROM: Level of nominal GDP (FRED :GDP series)
        # TO:   The cyclical component of log real GDP

        levels[:rGDP] = nominal_to_real(:GDP, levels)
        gdp = oneqtrpctchange(levels[:rGDP])
        _, levels[:GDP_cyc_comp] = hpfilter(gdp, 1600)
        levels[:GDP_cyc_comp]
        #=levels[:lGDP] = log.(levels[:rGDP])
        _, levels[:lGDP_cyc_comp] = hpfilter(levels[:lGDP], 1600)

        levels[:lGDP_cyc_comp]=#
    end

    gdp_rev_transform = logleveltopct_annualized

    observables[:obs_gdp] = Observable(:obs_gdp, [:GDP__FRED, population_mnemonic, :GDPDEF__FRED],
                                       gdp_fwd_transform, gdp_rev_transform,
                                       "Log Real GDP", "")

    ############################################################################
    ## 2. Core PCE Inflation
    ############################################################################

    pce_fwd_transform = function (levels)
        # FROM: Core PCE index
        # INTO: Approximate quarter-to-quarter percent change of Core PCE,
        # i.e. quarterly core pce inflation

        oneqtrpctchange(levels[:PCEPILFE])
    end

    pce_rev_transform = loggrowthtopct_annualized

    observables[:obs_corepce] = Observable(:obs_corepce, [:PCEPILFE__FRED],
                                           pce_fwd_transform, pce_rev_transform,
                                           "Core PCE Inflation",
                                           "Core PCE Inflation")

    ############################################################################
    ## 3. Nominal short-term interest rate (3 months)
    ############################################################################

    nominalrate_fwd_transform = function (levels)
        # FROM: Nominal effective federal funds rate (aggregate daily data at a
        #       quarterly frequency at an annual rate)
        # TO:   Nominal effective fed funds rate, at a quarterly rate

        annualtoquarter(levels[:DFF])
    end

    nominalrate_rev_transform = quartertoannual

    observables[:obs_nominalrate] = Observable(:obs_nominalrate, [:DFF__FRED],
                                               nominalrate_fwd_transform, nominalrate_rev_transform,
                                               "Nominal FFR",
                                               "Nominal Effective Fed Funds Rate")

    m.observable_mappings = observables
end
