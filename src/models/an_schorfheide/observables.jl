function init_observable_mappings!(m::AnSchorfheide)

    observables = OrderedDict{Symbol,Observable}()
    population_mnemonic = get(get_setting(m, :population_mnemonic))

    ############################################################################
    ## 1. Real GDP Growth
    ############################################################################
    gdp_fwd_transform = function (levels)
        # FROM: Level of GDP (from FRED)
        # TO: Quarter-to-quarter percent change of real GDP per capita

        levels[:temp] = percapita(m, :GDP, levels)
        gdp = 1000 * nominal_to_real(:temp, levels)
        hpadjust(oneqtrpctchange(gdp), levels)
    end

    gdp_rev_transform = logtopct_annualized_percapita

    observables[:obs_gdp] = Observable(:obs_gdp, [:GDP__FRED, population_mnemonic, :GDPCTPI__FRED],
                                       gdp_fwd_transform, gdp_rev_transform,
                                       "Real GDP Growth", "Real GDP Growth Per Capita")

    ############################################################################
    ## 2. CPI Inflation
    ############################################################################

    infl_fwd_transform = function (levels)
        # FROM: CPI urban consumers index (from FRED)
        # TO: Annualized quarter-to-quarter percent change of CPI index

        quartertoannual(oneqtrpctchange(levels[:CPIAUCSL]))
    end

    infl_rev_transform = logtopct_annualized

    observables[:obs_infl] = Observable(:obs_infl, [:CPIAUCSL__FRED],
                                        infl_fwd_transform, infl_rev_transform,
                                        "CPI Inflation",
                                        "CPI Inflation")

    ############################################################################
    ## 3. Nominal short-term interest rate (3 months)
    ############################################################################

    nominalrate_fwd_transform = function (levels)
        # FROM: Nominal effective federal funds rate (aggregate daily data at a
        #       quarterly frequency at an annual rate)
        # TO:   Nominal effective fed funds rate, at a quarterly rate

        levels[:DFF]
    end

    nominalrate_rev_transform = identity

    observables[:obs_nominalrate] = Observable(:obs_nominalrate, [:DFF__FRED],
                                               nominalrate_fwd_transform, nominalrate_rev_transform,
                                               "Nominal FFR",
                                               "Nominal Effective Fed Funds Rate")

    m.observable_mappings = observables
end
