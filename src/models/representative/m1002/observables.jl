function init_observable_mappings!(m::Model1002)

    observables = OrderedDict{Symbol,Observable}()
    population_mnemonic = get(get_setting(m, :population_mnemonic))

    ############################################################################
    ## 1. GDP
    ############################################################################
    gdp_fwd_transform =  function (levels)
        # FROM: Level of nominal GDP (FRED :GDP series)        # TO:   Quarter-to-quarter percent change of real, per-capita GDP, adjusted for population smoothing

        levels[!,:temp] = percapita(m, :GDP, levels)
        gdp = 1000 * nominal_to_real(:temp, levels)
        oneqtrpctchange(gdp)
    end

    gdp_rev_transform = loggrowthtopct_annualized_percapita

    ############################################################################
    ## 2. Hours per-capita
    ############################################################################

    hrs_fwd_transform =  function (levels)
        # FROM: Average weekly hours (AWHNONAG) & civilian employment (CE16OV)
        # TO:   log (3 * per-capita weekly hours / 100)
        # Note: Not sure why the 3 is there.
        # As of 11/6/2019: It is unclear where the 3 came from.
        #                  It is suspected that b/c the hours data is monthly,
        #                  the data had been processed so that you got the amount of hours
        #                  worked per month (on average), and multiplied by 3 to get a
        #                  quarterly guess of total hours worked.
        #                  Regardless, it doesn't matter for results b/c
        #                  it just shifts the value of Lmean, which doesn't appear anywhere
        #                  else and is just a normalizing constant for the model.
        #                  For consistency across the codebase, we maintain it as 3,
        #                  although it may be more accurate to do weeklyhours * 13.

        levels[!,:temp] = levels[!,:AWHNONAG] .* levels[!,:CE16OV]
        weeklyhours = percapita(m, :temp, levels)
        100*log.(3 * weeklyhours / 100)
    end

    hrs_rev_transform = logleveltopct_annualized_percapita

    if haskey(m.settings, :hours_first_observable)
        if get_setting(m, :hours_first_observable)
            observables[:obs_hours] = Observable(:obs_hours, [:AWHNONAG__FRED, :CE16OV__FRED],
                                                 hrs_fwd_transform, hrs_rev_transform,
                                                 "Hours Per Capita", "Log Hours Per Capita")
            observables[:obs_gdp] = Observable(:obs_gdp, [:GDP__FRED, population_mnemonic, :GDPDEF__FRED],
                                               gdp_fwd_transform, gdp_rev_transform,
                                               "Real GDP Growth", "Real GDP Growth Per Capita")
        else
            observables[:obs_gdp] = Observable(:obs_gdp, [:GDP__FRED, population_mnemonic, :GDPDEF__FRED],
                                               gdp_fwd_transform, gdp_rev_transform,
                                               "Real GDP Growth", "Real GDP Growth Per Capita")
            observables[:obs_hours] = Observable(:obs_hours, [:AWHNONAG__FRED, :CE16OV__FRED],
                                                 hrs_fwd_transform, hrs_rev_transform,
                                                 "Hours Per Capita", "Log Hours Per Capita")
        end
    else
        observables[:obs_gdp] = Observable(:obs_gdp, [:GDP__FRED, population_mnemonic, :GDPDEF__FRED],
                                           gdp_fwd_transform, gdp_rev_transform,
                                           "Real GDP Growth", "Real GDP Growth Per Capita")
        observables[:obs_hours] = Observable(:obs_hours, [:AWHNONAG__FRED, :CE16OV__FRED],
                                             hrs_fwd_transform, hrs_rev_transform,
                                             "Hours Per Capita", "Log Hours Per Capita")
    end

    ############################################################################
    ## 3. Wages
    ############################################################################
    if subspec(m) in ["ss16", "ss17"]
        laborshare_fwd_transform = function (levels)
            # FROM: Level of nominal GDP (FRED :GDP series), nominal compensation per hour,
            #       and average weekly hours (AWHNONAG) & civilian employment (CE16OV)
            # TO: log labor share of output, log(w * L / y)

            # Compute series 1: 100 * log(real wage) + log(hours) - log(real gdp)
            levels[!,:gdptemp] = percapita(m, :GDP, levels)
            gdp = 1000 * nominal_to_real(:gdptemp, levels)
            levels[!,:hourstemp] = levels[!,:AWHNONAG] .* levels[!,:CE16OV]
            weeklyhours = percapita(m, :hourstemp, levels)
            series1 = 100 * log.(nominal_to_real(:COMPNFB, levels)) +
                100 * log.(3 * weeklyhours / 100) - 100 * log.(gdp)

            # Compute base value and deviation of indices from initial value
            # log(BEA/NIPA nominal employee compensation / gdp) - log(series 1 in T0)
            laborshare_start_date = DSGE.firstdayofquarter(date_presample_start(m) -
                                                           Dates.Month(6))
            laborshare_end_date = date_mainsample_end(m)
            base_per = findfirst(get_setting(m, :laborshare_base_period) .==
                                 DSGE.get_quarter_ends(laborshare_start_date,
                                                       laborshare_end_date))
            series1_base = series1[base_per] # deviation of index from this base value
            series2_base = 100 * log(levels[!,:COE][base_per] /
                                     levels[!,:GDP][base_per]) # base labor share

            series1 .- series1_base .+ series2_base
        end

        laborshare_rev_transform = logleveltopct_annualized_percapita

        observables[:obs_laborshare] = Observable(:obs_laborshare, [:COMPNFB__FRED, :GDPDEF__FRED,
                                                                   :COE__FRED],
                                             laborshare_fwd_transform, laborshare_rev_transform,
                                             "Log Labor Share",
                                             "Quarterly Log Labor Share of Real GDP")
    else
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
    end
    ############################################################################
    ## 4. GDP Deflator
    ############################################################################

    gdpdeflator_fwd_transform =  function (levels)
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
    ## 5. Core PCE Inflation
    ############################################################################

    pce_fwd_transform = function (levels)
        # FROM: Core PCE index
        # INTO: Approximate quarter-to-quarter percent change of Core PCE,
        # i.e. quarterly core pce inflation

        oneqtrpctchange(levels[!,:PCEPILFE])
    end

    pce_rev_transform = loggrowthtopct_annualized

    observables[:obs_corepce] = Observable(:obs_corepce, [:PCEPILFE__FRED],
                                           pce_fwd_transform, pce_rev_transform,
                                           "Core PCE Inflation",
                                           "Core PCE Inflation")


    ############################################################################
    ## 6. Nominal short-term interest rate (3 months)
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
    ## 7. Consumption
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
    ## 8. Investment growth per capita
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
    ## 9. Spread: BAA-10yr TBill
    ############################################################################

    spread_fwd_transform = function (levels)
        # FROM: Baa corporate bond yield (percent annualized), and 10-year
        #       treasury note yield (percent annualized)
        # TO:   Baa yield - 10T yield spread at a quarterly rate
        # Note: Moody's corporate bond yields on the H15 are based on corporate
        #       bonds with remaining maturities of at least 20 years.
        #       The Moody's series (BAA) ends partway through 2016-Q4. Hence
        #       beginning in 2016-Q4, we replace it with a similar series from
        #       Bank of America (BAMLC8A0C15PYEY).

        splice_date   = quartertodate("2016-Q4") # quarter at which we start using new series
        old_series    = levels[levels[!,:date] .<  splice_date, :BAA]
        new_series    = levels[levels[!,:date] .>= splice_date, :BAMLC8A0C15PYEY]
        levels[!,:temp] = vcat(old_series, new_series)

        annualtoquarter(levels[!,:temp] - levels[!,:GS10])
    end

    spread_rev_transform = quartertoannual

    observables[:obs_spread] = Observable(:obs_spread, [:BAA__FRED, :BAMLC8A0C15PYEY__FRED, :GS10__FRED],
                                          spread_fwd_transform, spread_rev_transform,
                                          "BAA - 10yr Treasury Spread",
                                          "BAA - 10yr Treasury Spread")


    ############################################################################
    # 10. Long term inflation expectations
    ############################################################################

    longinflation_fwd_transform = function (levels)
        # FROM: SPF: 10-Year average yr/yr CPI inflation expectations (annual percent)
        # TO:   FROM, less 0.5
        # Note: We subtract 0.5 because 0.5% inflation corresponds to
        #       the assumed long-term rate of 2 percent inflation, but the
        #       data are measuring expectations of actual inflation.

        annualtoquarter(levels[!,:ASACX10]  .- 0.5)
    end

    longinflation_rev_transform = loggrowthtopct_annualized

    observables[:obs_longinflation] = Observable(:obs_longinflation, [:ASACX10__DLX],
                                                 longinflation_fwd_transform, longinflation_rev_transform,
                                                 "Long term inflation expectations",
                                                 "10-year average yr/yr CPI inflation expectations")


    ############################################################################
    # 11. Long rate (10-year, zero-coupon)
    ############################################################################
    longrate_fwd_transform = function (levels)
        # FROM: pre-computed long rate at an annual rate
        # TO:   10T yield at a quarterly rate

        annualtoquarter(levels[!,:FYCCZA])
    end

    longrate_rev_transform = quartertoannual

    observables[:obs_longrate] = Observable(:obs_longrate, [:FYCCZA__DLX],
                                            longrate_fwd_transform, longrate_rev_transform,
                                            "Long term interest rate expectations",
                                            "10T yield")


    ############################################################################
    # 12. Fernald TFP
    ############################################################################

    tfp_rev_transform = quartertoannual

    if subspec(m) in ["ss15", "ss16"]
        tfp_fwd_transform =  function (levels)
            # FROM: Fernald's adjusted TFP series
            # TO:   De-meaned adjusted TFP series, adjusted by Fernald's estimated alpha
            # Note: We only want to calculate the mean of unadjusted/adjusted TFP over the
            #       periods between date_presample_start(m) - 1 and
            #       date_mainsample_end(m), though we may end up transforming
            #       additional periods of data.

            start_date = Dates.lastdayofquarter(date_presample_start(m) - Dates.Month(3))
            end_date   = date_mainsample_end(m)
            date_range = start_date .<= levels[1:end, :date] .<= end_date
            tfp_unadj_inrange = levels[date_range, :TFPMQ]

            tfp_unadj      = levels[!,:TFPMQ]
            tfp_unadj_inrange_nonmissing = tfp_unadj_inrange[.!ismissing.(tfp_unadj_inrange)]
            tfp_unadj_inrange_nonmissing = tfp_unadj_inrange_nonmissing[.!isnan.(tfp_unadj_inrange_nonmissing)]
            tfp_unadj_mean = isempty(tfp_unadj_inrange_nonmissing) ? missing : mean(tfp_unadj_inrange_nonmissing)
            (tfp_unadj .- tfp_unadj_mean) ./ (4*(1 .- levels[!,:TFPJQ]))
        end

        observables[:obs_tfp] = Observable(:obs_tfp, [:TFPMQ__tfputiladj, :TFPJQ__tfputiladj],
                                           tfp_fwd_transform, tfp_rev_transform,
                                           "Total Factor Productivity Growth (Fernald)",
                                           "Fernald's TFP, adjusted by Fernald's estimated alpha and utilization capacity")
    else
        tfp_fwd_transform =  function (levels)
            # FROM: Fernald's unadjusted TFP series
            # TO:   De-meaned unadjusted TFP series, adjusted by Fernald's estimated alpha
            # Note: We only want to calculate the mean of unadjusted/adjusted TFP over the
            #       periods between date_presample_start(m) - 1 and
            #       date_mainsample_end(m), though we may end up transforming
            #       additional periods of data.

            start_date = Dates.lastdayofquarter(date_presample_start(m) - Dates.Month(3))
            end_date   = date_mainsample_end(m)
            date_range = start_date .<= levels[1:end, :date] .<= end_date
            tfp_unadj_inrange = levels[date_range, :TFPKQ]

            tfp_unadj      = levels[!,:TFPKQ]
            tfp_unadj_inrange_nonmissing = tfp_unadj_inrange[.!ismissing.(tfp_unadj_inrange)]
            tfp_unadj_inrange_nonmissing = tfp_unadj_inrange_nonmissing[.!isnan.(tfp_unadj_inrange_nonmissing)]
            tfp_unadj_mean = isempty(tfp_unadj_inrange_nonmissing) ? missing : mean(tfp_unadj_inrange_nonmissing)
            (tfp_unadj .- tfp_unadj_mean) ./ (4*(1 .- levels[!,:TFPJQ]))
        end

        observables[:obs_tfp] = Observable(:obs_tfp, [:TFPKQ__DLX, :TFPJQ__DLX],
                                           tfp_fwd_transform, tfp_rev_transform,
                                           "Total Factor Productivity Growth (Fernald)",
                                           "Fernald's TFP, adjusted by Fernald's estimated alpha")
    end

    ############################################################################
    # 13. GDI
    ############################################################################
    gdi_fwd_transform = function (levels)
        # FROM: level of real per-capita income (FRED mnemonic: GDI)
        # TO:   approximate quarter-to-quarter percent change of real, per-capita income, adjusted
        #       for population smoothing

        levels[!,:temp] = percapita(m, :GDI, levels)
        gdi = 1000 * nominal_to_real(:temp, levels)
        oneqtrpctchange(gdi)
    end

    gdi_rev_transform = loggrowthtopct_annualized_percapita

    observables[:obs_gdi] = Observable(:obs_gdi, [:GDI__FRED],
                                       gdi_fwd_transform, gdi_rev_transform,
                                       "Real GDI Growth",
                                       "Real GDI Growth Per Capita")

    ############################################################################
    # Columns 14 - 14 + n_mon_anticipated_shocks
    ############################################################################

    for i = 1:n_mon_anticipated_shocks(m)
        # FROM: OIS expectations of $i-period-ahead interest rates at a quarterly rate
        # TO:   Same

        ant_fwd_transform = function (levels)
            levels[:, Symbol("ant$i")]
        end

        ant_rev_transform = quartertoannual

        observables[Symbol("obs_nominalrate$i")] = Observable(Symbol("obs_ant$i"), [Symbol("ant$(i)__OIS")],
                                                      ant_fwd_transform, ant_rev_transform,
                                                      "Anticipated Shock $i",
                                                      "$i-period ahead anticipated monetary policy shock")
    end

    ############################################################################
    # Other anticipated data
    ############################################################################
    if haskey(get_settings(m), :add_anticipated_obs_gdp)
        if get_setting(m, :add_anticipated_obs_gdp)
            antgdp_name = haskey(get_settings(m), :filename_anticipated_obs_gdp) ?
                get_setting(m, :filename_anticipated_obs_gdp) : "ANTGDP"
            for i = 1:get_setting(m, :n_anticipated_obs_gdp)
                # FROM: Some source for expectations off i-period ahead GDP growth
                # TO:   Same

                ant_fwd_transform = function (levels)
                    levels[:, Symbol("antgdp$i")]
                end

                ant_rev_transform = loggrowthtopct_annualized_percapita

                observables[Symbol("obs_gdp$i")] = Observable(Symbol("obs_antgdp$i"), [Symbol("antgdp$(i)__$(antgdp_name)")],
                                                                      ant_fwd_transform, ant_rev_transform,
                                                                      "Anticipated GDP Growth $i",
                                                                      "$i-period ahead anticipated GDP growth")
            end
        end
    end

    ############################################################################
    # Pseudo-data to implement Flexible AIT
    ############################################################################
    if (haskey(m.settings, :add_initialize_pgap_ygap_pseudoobs) ? get_setting(m, :add_initialize_pgap_ygap_pseudoobs) : false)
        pgap_fwd_transform = function (levels)
            levels[:, :pgap]
        end


        observables[:obs_pgap] = Observable(:obs_pgap, [:pgap__INITFLEXAIT],
                                            pgap_fwd_transform, identity,
                                            "Average Inflation Gap",
                                            "Average Inflation Gap from Target")

        ygap_fwd_transform = function (levels)
            levels[:, :ygap]
        end

        observables[:obs_ygap] = Observable(:obs_ygap, [:ygap__INITFLEXAIT],
                                            ygap_fwd_transform, identity,
                                            "Average Output Gap",
                                            "Average Output Gap from Target")
    end

    if haskey(m.settings, :first_observable)
        new_observables = OrderedDict{Symbol,Observable}()
        first_obs = get_setting(m, :first_observable)
        new_observables[first_obs] = observables[first_obs]
        for (k,v) in observables
            if k != first_obs
                new_observables[k] = v
            end
        end
        observables = new_observables
    end
    if haskey(m.settings, :last_observable)
        new_observables = OrderedDict{Symbol,Observable}()
        last_obs = get_setting(m, :last_observable)
        for (k,v) in observables
            if k != last_obs
                new_observables[k] = v
            end
        end
        new_observables[last_obs] = observables[last_obs]
        observables = new_observables
    end

    # Needed to implement measurement equation correctly
    m <= Setting(:forward_looking_observables,
                 vcat([:obs_longinflation, :obs_longrate],
                      [Symbol("obs_nominalrate$i") for i in 1:n_mon_anticipated_shocks(m)],
                      haskey(get_settings(m), :add_anticipated_obs_gdp) && get_setting(m, :add_anticipated_obs_gdp) ?
                      [Symbol("obs_gdp$i") for i in 1:get_setting(m, :n_anticipated_obs_gdp)] : []))

    m.observable_mappings = observables
end
