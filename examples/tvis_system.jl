using DSGE, ModelConstructors, Dates, Plots
fn = dirname(@__FILE__)

# This script shows how to calculate state space systems
# with time-varying information sets and shows how they
# affect forecasts. See the documentation
# for mathematical and theoretical clarifications on
# what we mean precisely by time-varying information sets.

# What do you want to do?
calculate_system = true
modal_forecast   = true
create_plots     = true

## Set up

# Initialize model w/desired dates
m = Model1002("ss10")
usual_model_settings!(m, "200930"; cdvt = "200930", fcast_date = quartertodate("2020-Q3"))
m <= Setting(:use_population_forecast, false)
m <= Setting(:saveroot, joinpath("$(fn)", "..", "save"))
m <= Setting(:dataroot, joinpath("$(fn)", "..", "save", "input_data"))
m <= Setting(:date_forecast_end, quartertodate("2029-Q4"))

if calculate_system
    # Set up two different sets of equilibrium conditions,
    # one of which corresponds to an expected temporary ZLB
    # from 2020:Q3 to 2023:Q4 and another which
    # corresponds to an expected temporary ZLB from
    # 2020:Q3 to 2024:Q4 (one additional year).
    m <= Setting(:replace_eqcond, true) # turn on replace_eqcond and gensys2
    m <= Setting(:gensys2, true)        # to allow temporary alternative policies (namely the ZLB)
    regime_eqcond_info_shortzlb = Dict{Int, EqcondEntry}()
    replace_eqcond_func_dict_shortzlb = Dict()
    for i in 2:(length(DSGE.quarter_range(Date(2020, 9, 30), Date(2023, 12, 31))) + 1)
        # replace_eqcond_func_dict_shortzlb[i] = DSGE.zero_rate_replace_eq_entries
        regime_eqcond_info_shortzlb[i] = EqcondEntry(zero_rate())
    end
    regime_eqcond_info_longzlb = Dict{Int, EqcondEntry}()
    replace_eqcond_func_dict_longzlb = Dict()
    for i in 2:(length(DSGE.quarter_range(Date(2020, 9, 30), Date(2024, 12, 31))) + 1)
        # replace_eqcond_func_dict_longzlb[i] = DSGE.zero_rate_replace_eq_entries
        regime_eqcond_info_longzlb[i] = EqcondEntry(zero_rate())
    end
    m <= Setting(:zero_rate_zlb_value, 0.)
    m <= Setting(:temporary_altpolicy_names, [:zero_rate])

    # Set up regime-switching dates and indices
    regime_dates = Dict{Int, Date}(1 => date_presample_start(m))
    for (i, d) in enumerate(DSGE.quarter_range(Date(2020, 9, 30), Date(2025, 3, 31)))
        # Add regimes for each quarter up to 2025:Q1 b/c the same :regime_dates dictionary
        # will be used when calculating the state-space system for the shorter ZLB policy
        # and the longer ZLB policy. If we added regimes only up to 2024:Q1, then
        # the longer ZLB policy will run into an error.
        regime_dates[i + 1] = d
    end
    m <= Setting(:regime_dates, regime_dates)
    m <= Setting(:regime_switching, true)
    setup_regime_switching_inds!(m)

    ## Calculate state-space system with time-varying information sets

    # First, let's calculate an example where people learn about the short ZLB
    # starting in 2020:Q3. In 2020:Q2 (and before), people do NOT expect a ZLB.
    reg_2024Q1 = get_setting(m, :n_regimes) - 4 # This is the regime corresponding to 2024:Q1
    @assert get_setting(m, :regime_dates)[reg_2024Q1] == Date(2024, 3, 31)

    # How the time-varying information set is constructed:
    # In regime 1, agents only need to use the regime 1 transition equation
    # From 2020:Q3 to 2024:Q1, people need to use the transition equations for 2020:Q3 to 2024:Q1
    # to construct the measurement equation. In 2024:Q2 onward, people again only need to use
    # the current regime's transition equation. Alternatively, since the transition equation
    # in 2024:Q1 is the same as the one in 2024:Q2 - 2025:Q1, we could have also assumed that,
    # during 2020:Q3 to 2024:Q1, people use the transition equation from 2020:Q3 to 2025:Q1.
    # This alternative approach is mathematically equivalent.
    m <= Setting(:tvis_information_set,
                 vcat([1:1],
                      [i:reg_2024Q1 for i in 2:reg_2024Q1],
                      [[i:i] for i in (reg_2024Q1 + 1):get_setting(m, :n_regimes)]...))
    # m <= Setting(:replace_eqcond_func_dict, replace_eqcond_func_dict_shortzlb)
    m <= Setting(:regime_eqcond_info, regime_eqcond_info_shortzlb)
    system_shortzlb = compute_system(m; tvis = true)

    # We could use also these settings to generate the same results.
    # m <= Setting(:tvis_replace_eqcond_func_dict, [replace_eqcond_func_dict_shortzlb])
    # m <= Setting(:tvis_select_system, ones(Int, 7))

    # Note that if :tvis_replace_eqcond_func_dict is already a setting in m,
    # then we use this setting to determine the dictionary of functions for
    # :replace_eqcond_func_dict rather than what is already stored
    # in the setting :replace_eqcond_func_dict.

    # Now, let's calculate an example where people learn about the short ZLB
    # starting in 2020:Q3, but in 2021:Q2, people learn that the long ZLB
    # will go into action.

    # First, we add info about the equilibrium conditions
    reg_2021Q1 = DSGE.subtract_quarters(Date(2021, 3, 31), Date(2020, 9, 30)) + 2
    # m <= Setting(:tvis_replace_eqcond_func_dict, [replace_eqcond_func_dict_shortzlb, replace_eqcond_func_dict_longzlb])
    m <= Setting(:tvis_regime_eqcond_info, [regime_eqcond_info_shortzlb, regime_eqcond_info_longzlb])
    m <= Setting(:tvis_select_system, vcat(ones(Int, reg_2021Q1), # Use set 1 of conditions for regimes up to (and incl.) 2021:Q1
                                           fill(2, get_setting(m, :n_regimes) - reg_2021Q1))) # Use set 2 for remaining regimes

    # Second, we add info about the information sets
    m <= Setting(:tvis_information_set,
                 vcat([1:1], # Regime 1 doesn't change
                      [i:reg_2024Q1 for i in 2:reg_2021Q1], # This is the same as before, but the info set changes after 2021:Q1
                      [i:get_setting(m, :n_regimes) for i in (reg_2021Q1 + 1):get_setting(m, :n_regimes)])) # Now we use regimes to 2025:Q1 for the info set
    system_shortlong_zlb = compute_system(m; tvis = true)

    # Some tests that should pass
    for i in 1:reg_2021Q1
        @assert system_shortzlb[i, :ZZ] == system_shortlong_zlb[i, :ZZ]
        @assert system_shortzlb[i, :TTT] == system_shortlong_zlb[i, :TTT]
        @assert system_shortzlb[i, :ZZ] != system_shortlong_zlb[(reg_2021Q1 + 1), :ZZ]
        @assert system_shortzlb[i, :TTT] != system_shortlong_zlb[(reg_2021Q1 + 1), :TTT]
    end
    for i in (reg_2021Q1 + 1):(get_setting(m, :n_regimes) - 2)
        @assert system_shortzlb[i, :ZZ] != system_shortlong_zlb[i, :ZZ]
        @assert system_shortzlb[i, :TTT] != system_shortlong_zlb[i, :TTT]
    end
    @assert system_shortzlb[get_setting(m, :n_regimes) - 1, :ZZ] ≈ # b/c transition equation matters only for
        system_shortlong_zlb[get_setting(m, :n_regimes) - 1, :ZZ]  # forward expectations in measurement equation
    @assert system_shortzlb[get_setting(m, :n_regimes), :ZZ] ==
        system_shortlong_zlb[get_setting(m, :n_regimes), :ZZ]
    @assert system_shortzlb[get_setting(m, :n_regimes), :TTT] ==
        system_shortlong_zlb[get_setting(m, :n_regimes), :TTT]
end

if modal_forecast
    # Set up in case calculate_system is not true
    if !calculate_system
        m <= Setting(:replace_eqcond, true)
        m <= Setting(:gensys2, true)
        regime_eqcond_info_shortzlb = Dict{Int, EqcondEntry}()
        replace_eqcond_func_dict_shortzlb = Dict()
        for i in 2:(length(DSGE.quarter_range(Date(2020, 9, 30), Date(2023, 12, 31))) + 1)
            # replace_eqcond_func_dict_shortzlb[i] = DSGE.zero_rate_replace_eq_entries
            regime_eqcond_info_shortzlb[i] = EqcondEntry(zero_rate())
        end
        regime_eqcond_info_longzlb = Dict{Int, EqcondEntry}()
        replace_eqcond_func_dict_longzlb = Dict()
        for i in 2:(length(DSGE.quarter_range(Date(2020, 9, 30), Date(2024, 12, 31))) + 1)
            # replace_eqcond_func_dict_longzlb[i] = DSGE.zero_rate_replace_eq_entries
            regime_eqcond_info_longzlb[i] = EqcondEntry(zero_rate())
        end

        regime_dates = Dict{Int, Date}(1 => date_presample_start(m))
        for (i, d) in enumerate(DSGE.quarter_range(Date(2020, 9, 30), Date(2025, 3, 31)))
            regime_dates[i + 1] = d
        end
        m <= Setting(:zero_rate_zlb_value, 0.)
        m <= Setting(:temporary_altpolicy_names, [:zero_rate])

        m <= Setting(:regime_dates, regime_dates)
        m <= Setting(:regime_switching, true)
        setup_regime_switching_inds!(m)
    end

    # Set up two modal forecasts, to show the effect of extending the ZLB
    reg_2024Q1 = get_setting(m, :n_regimes) - 4
    m <= Setting(:tvis_information_set, vcat([1:1], [i:reg_2024Q1 for i in 2:reg_2024Q1],
                                             [[i:i] for i in (reg_2024Q1 + 1):get_setting(m, :n_regimes)]...))
    # m <= Setting(:tvis_replace_eqcond_func_dict, [replace_eqcond_func_dict_shortzlb])
    m <= Setting(:tvis_regime_eqcond_info, [regime_eqcond_info_shortzlb])
    # m <= Setting(:replace_eqcond_func_dict, replace_eqcond_func_dict_shortzlb)
    m <= Setting(:regime_eqcond_info, regime_eqcond_info_shortzlb)

    usual_model_forecast(m, :mode, :none, [:histobs, :forecastobs, :histpseudo, :forecastpseudo];
                         forecast_string = "short", params = [x.value for x in m.parameters],
                         check_empty_columns = false)
    mbobs_short = read_mb(m, :mode, :none, :forecastobs; forecast_string = "short")
    mbpse_short = read_mb(m, :mode, :none, :forecastpseudo; forecast_string = "short")

    reg_2021Q1 = DSGE.subtract_quarters(Date(2021, 3, 31), Date(2020, 9, 30)) + 2
    # m <= Setting(:tvis_replace_eqcond_func_dict, [replace_eqcond_func_dict_shortzlb, replace_eqcond_func_dict_longzlb])
    m <= Setting(:tvis_regime_eqcond_info, [regime_eqcond_info_shortzlb, regime_eqcond_info_longzlb])
    m <= Setting(:tvis_select_system, vcat(ones(Int, reg_2021Q1), fill(2, get_setting(m, :n_regimes) - reg_2021Q1)))
    m <= Setting(:tvis_information_set, vcat([1:1], [i:reg_2024Q1 for i in 2:reg_2021Q1],
                                             [i:get_setting(m, :n_regimes) for i in (reg_2021Q1 + 1):get_setting(m, :n_regimes)]))
    usual_model_forecast(m, :mode, :none, [:histobs, :forecastobs, :histpseudo, :forecastpseudo];
                         forecast_string = "shortlong", params = [x.value for x in m.parameters],
                         check_empty_columns = false)
    mbobs_shortlong = read_mb(m, :mode, :none, :forecastobs; forecast_string = "shortlong")
    mbpse_shortlong = read_mb(m, :mode, :none, :forecastpseudo; forecast_string = "shortlong")

    if create_plots
        plot_dicts = Dict()
        dates = [DSGE.iterate_quarters(d, -1) for d in mbobs_short.means[!, :date]]

        plot_dicts[:obs_gdp] = plot()
        plot!(dates, mbobs_shortlong.means[!, :obs_gdp], color = :black, label = "2024:Q4 ZLB", linewidth = 3)
        plot!(dates, mbobs_short.means[!, :obs_gdp], color = :blue, label = "2023:Q4 ZLB", linewidth = 3)

        plot_dicts[:obs_corepce] = plot()
        plot!(dates, mbobs_shortlong.means[!, :obs_corepce], color = :black, label = "2024:Q4 ZLB", linewidth = 3)
        plot!(dates, mbobs_short.means[!, :obs_corepce], color = :blue, label = "2023:Q4 ZLB", linewidth = 3)

        plot_dicts[:obs_nominalrate] = plot()
        plot!(dates, mbobs_shortlong.means[!, :obs_nominalrate], color = :black, label = "2024:Q4 ZLB", linewidth = 3)
        plot!(dates, mbobs_short.means[!, :obs_nominalrate], color = :blue, label = "2023:Q4 ZLB", linewidth = 3)

        plot_dicts[:obs_longinflation] = plot()
        plot!(dates, mbobs_shortlong.means[!, :obs_longinflation], color = :black, label = "2024:Q4 ZLB", linewidth = 3)
        plot!(dates, mbobs_short.means[!, :obs_longinflation], color = :blue, label = "2023:Q4 ZLB", linewidth = 3)

        plot_dicts[:obs_spread] = plot()
        plot!(dates, mbobs_shortlong.means[!, :obs_spread], color = :black, label = "2024:Q4 ZLB", linewidth = 3)
        plot!(dates, mbobs_short.means[!, :obs_spread], color = :blue, label = "2023:Q4 ZLB", linewidth = 3)

        plot_dicts[:NaturalRate] = plot()
        plot!(dates, mbpse_shortlong.means[!, :NaturalRate], color = :black, label = "2024:Q4 ZLB", linewidth = 3)
        plot!(dates, mbpse_short.means[!, :NaturalRate], color = :blue, label = "2023:Q4 ZLB", linewidth = 3)
    end
end
