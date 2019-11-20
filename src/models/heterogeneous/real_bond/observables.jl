function init_observable_mappings!(m::RealBond)

    observables = OrderedDict{Symbol,Observable}()
    population_mnemonic = get(get_setting(m, :population_mnemonic))

    ############################################################################
    ## 1. GDP
    ############################################################################
    gdp_fwd_transform = function (levels)
        # FROM: Level of nominal GDP (FRED :GDP series)
        # TO:   The cyclical component of log real GDP

        levels[:rGDP] = nominal_to_real(:GDP, levels)
        levels[:lGDP] = log.(levels[:rGDP])
        _, levels[:lGDP_cyc_comp] = hpfilter(levels[:lGDP], 1600)

        levels[:lGDP_cyc_comp]
    end

    gdp_rev_transform = logleveltopct_annualized

    observables[:obs_gdp] = Observable(:obs_gdp, [:GDP__FRED, population_mnemonic, :GDPDEF__FRED],
                                       gdp_fwd_transform, gdp_rev_transform,
                                       "Log Real GDP", "")

    m.observable_mappings = observables
end
