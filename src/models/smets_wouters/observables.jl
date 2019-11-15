function init_observable_mappings!(m::SmetsWouters)

    observables = OrderedDict{Symbol,Observable}()
    population_mnemonic = get(get_setting(m, :population_mnemonic))

    ############################################################################
    ## 1. GDP
    ############################################################################
    gdp_fwd_transform = function (levels)
        # FROM: Level of nominal GDP (FRED :GDP series)
        # TO:   Quarter-to-quarter percent change of real, per-capita GDP, adjusted for population smoothing

        levels[!,:temp] = percapita(m, :GDP, levels)
        gdp = 1000 * nominal_to_real(:temp, levels)
        oneqtrpctchange(gdp)
    end

    gdp_rev_transform = loggrowthtopct_annualized_percapita

    observables[:obs_gdp] = Observable(:obs_gdp, [:GDP__FRED, population_mnemonic, :GDPDEF__FRED],
                                       gdp_fwd_transform, gdp_rev_transform,
                                       "Real GDP Growth", "Real GDP Growth Per Capita")

    ############################################################################
    ## 2. Hours per-capita
    ############################################################################

    hrs_fwd_transform = function (levels)
        # FROM: Average weekly hours (AWHNONAG) & civilian employment (CE16OV)
        # TO:   log (3 * per-capita weekly hours / 100)
        # Note: Not sure why the 3 is there.

        levels[!,:temp] = levels[!,:AWHNONAG] .* levels[!,:CE16OV]
        weeklyhours = percapita(m, :temp, levels)
        100*log.(3 * weeklyhours / 100)
    end

    hrs_rev_transform = logleveltopct_annualized_percapita

    observables[:obs_hours] = Observable(:obs_hours, [:AWHNONAG__FRED, :CE16OV__FRED],
                                         hrs_fwd_transform, hrs_rev_transform,
                                         "Hours Per Capita", "Log Hours Per Capita")

    ############################################################################
    ## 3. Wages
    ############################################################################

    wages_fwd_transform = function (levels)
        # FROM: Nominal compensation per hour (:COMPNFB from FRED)
        # TO: quarter to quarter percent change of real compensation (using GDP deflator)

        oneqtrpctchange(nominal_to_real(:COMPNFB, levels))
    end

    wages_rev_transform = loggrowthtopct_annualized

    observables[:obs_wages] = Observable(:obs_wages, [:COMPNFB__FRED, :GDPDEF__FRED],
                                         wages_fwd_transform, wages_rev_transform,
                                         "Percent Change in Wages",
                                         "Q-to-Q Percent Change of Real Compensation (using GDP deflator)")

    ############################################################################
    ## 4. GDP Deflator
    ############################################################################

    gdpdeflator_fwd_transform = function (levels)
        # FROM: GDP deflator (index)
        # TO:   Approximate quarter-to-quarter percent change of gdp deflator,
        #       i.e.  quarterly gdp deflator inflation

        oneqtrpctchange(levels[!,:GDPDEF])
    end


    gdpdeflator_rev_transform = loggrowthtopct_annualized

    observables[:obs_gdpdeflator] = Observable(:obs_gdpdeflator, [:GDPDEF__FRED],
                                               gdpdeflator_fwd_transform, gdpdeflator_rev_transform,
                                               "GDP Deflator",
                                               "Q-to-Q Percent Change of GDP Deflator")

    ############################################################################
    ## 5. Nominal short-term interest rate (3 months)
    ############################################################################

    nominalrate_fwd_transform = function (levels)
        # FROM: Nominal effective federal funds rate (aggregate daily data at a
        #       quarterly frequency at an annual rate)
        # TO:   Nominal effective fed funds rate, at a quarterly rate

        annualtoquarter(levels[!,:DFF])
    end

    nominalrate_rev_transform = quartertoannual

    observables[:obs_nominalrate] = Observable(:obs_nominalrate, [:DFF__FRED],
                                               nominalrate_fwd_transform, nominalrate_rev_transform,
                                               "Nominal FFR",
                                               "Nominal Effective Fed Funds Rate")

    ############################################################################
    ## 6. Consumption
    ############################################################################

    consumption_fwd_transform = function (levels)
        # FROM: Nominal consumption
        # TO:   Real consumption, approximate quarter-to-quarter percent change,
        #       per capita, adjusted for population filtering

        levels[!,:temp] = percapita(m, :PCE, levels)
        cons = 1000 * nominal_to_real(:temp, levels)
        oneqtrpctchange(cons)
    end

    consumption_rev_transform = loggrowthtopct_annualized_percapita

    observables[:obs_consumption] = Observable(:obs_consumption, [:PCE__FRED, population_mnemonic],
                                               consumption_fwd_transform, consumption_rev_transform,
                                               "Consumption growth per capita",
                                               "Consumption growth adjusted for population filtering")

    ############################################################################
    ## 7. Investment growth per capita
    ############################################################################

    investment_fwd_transform = function (levels)
        # FROM: Nominal investment
        # INTO: Real investment, approximate quarter-to-quarter percent change,
        #       per capita, adjusted for population filtering

        levels[!,:temp] = percapita(m, :FPI, levels)
        inv = 10000 * nominal_to_real(:temp, levels)
        oneqtrpctchange(inv)
    end

    investment_rev_transform  = loggrowthtopct_annualized_percapita

    observables[:obs_investment] = Observable(:obs_investment, [:FPI__FRED, population_mnemonic],
                                              investment_fwd_transform, investment_rev_transform,
                                              "Real Investment per capita",
                                              "Real investment per capita, adjusted for population filtering")

    ############################################################################
    # 8 - (8 + n_anticipated_shocks)
    ############################################################################

    for i = 1:n_anticipated_shocks(m)
        # FROM: OIS expectations of $i-period-ahead interest rates at a quarterly rate
        # TO:   Same

        ant_fwd_transform = function (levels)
            levels[!, Symbol("ant$i")]
        end

        ant_rev_transform = quartertoannual

        observables[Symbol("obs_nominalrate$i")] = Observable(Symbol("obs_ant$i"), [Symbol("ant$(i)__OIS")],
                                                      ant_fwd_transform, ant_rev_transform,
                                                      "Anticipated Shock $i",
                                                      "$i-period ahead anticipated monetary policy shock")
    end

    m.observable_mappings = observables
end
