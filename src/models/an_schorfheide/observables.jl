function init_observable_mappings!(m::AnSchorfheide)

    observables = OrderedDict{Symbol,Observable}()

    ############################################################################
    ## 1. Real GDP Growth
    ############################################################################
    gdp_fwd_transform = function (levels)
        # FROM: Level of real GDP (from FRED)
        # TO: Quarter-to-quarter percent change of real GDP

        levels[:temp] = percapita(m, :GDP, levels)
        gdp = 1000 * nominal_to_real(:temp, levels)
        hpadjust(oneqtrpctchange(gdp), levels)
    end

    gdp_rev_transform = DSGE.logtopct_annualized_percapita

    observables[:obs_gdp] = Observable(:obs_gdp, [:GDP__FRED, :CNP16OV__FRED, :GDPCTPI__FRED],
                                       gdp_fwd_transform, gdp_rev_transform,
                                       "Real GDP Growth", "Real GDP Growth Per Capita")

    ############################################################################
    ## 2. Core PCE Inflation
    ############################################################################

    pce_fwd_transform = function (levels)
        # FROM: Core PCE index (from FRED)
        # TO: Quarter-to-quarter percent change of core PCE, i.e. quarterly inflation
        oneqtrpctchange(levels[:PCEPILFE])
    end

    pce_rev_transform = logtopct_annualized

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
