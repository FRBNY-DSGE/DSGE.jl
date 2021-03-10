using ClusterManagers, DSGE, ModelConstructors, FileIO, Plots, Dates, OrderedCollections, Distributed
using Plots.PlotMeasures
GR.inline("pdf")
fn = dirname(@__FILE__)

# This script shows how to use regime-switching
# to enact a temporary ZLB with imperfect awareness, followed by an
# alternative policy with imperfect awareness (in this case, AIT).
# We sometimes also call policies with imperfect awareness as "uncertain" policies.
# This script demonstrates a specific case of the general regime switching code
# demonstrated in regime_switching.jl.

# What do you want to do?
zlb_altpolicy     = true   # run a forecast
make_plots        = false  # Make plots
save_plots        = false  # Save plots to figurespath(m, "forecast")
add_workers       = false  # Run in parallel
n_workers         = 10

pol_str = "Uncertain ZLB, Uncertain AIT"

# Initialize model objects and desired settings
custom_settings = Dict{Symbol, Setting}(:n_mon_anticipated_shocks =>
                                        Setting(:n_mon_anticipated_shocks, 6, "Number of anticipated policy shocks"),
                                        :add_altpolicy_pgap => Setting(:add_altpolicy_pgap, true),
                                        :add_altpolicy_yagp => Setting(:add_altpolicy_ygap, true),
                                        :flextible_ait_policy_change =>
                                        Setting(:flexible_ait_policy_change, false))
m = Model1002("ss10"; custom_settings = custom_settings)   # We will directly construct the matrix of parameter draws
m <= Setting(:flexible_ait_2020Q3_policy_change, false)

data_vint = "201117"
cond_vint = "201119"
fcast_date = DSGE.quartertodate("2020-Q4")
cond_date = DSGE.quartertodate("2020-Q4")
date_fcast_end = iterate_quarters(fcast_date, 40)

m10 = Model1002("ss10") # b/c loading from a saved estimation file has not been fully implemented

for models in [m, m10]
    # Not strictly necessary to have the same settings, but it makes it easier for us
    usual_model_settings!(m, "201117"; cdvt = "201119", fcast_date = quartertodate("2020-Q4"))
    m <= Setting(:use_population_forecast, false)
    m <= Setting(:saveroot, joinpath("$(fn)", "..", "save"))
    m <= Setting(:dataroot, joinpath("$(fn)", "..", "save", "input_data"))
    m <= Setting(:date_forecast_end, date_fcast_end)
    m <= Setting(:forecast_block_size, 40)
    if get_setting(m, :sampling_method) == :SMC
        m <= Setting(:forecast_jstep, 1)
    end
end

# Set up parameters for regime-switching.
overrides = forecast_input_file_overrides(m10)
overrides[:full] = "$(fn)/../test/reference/mhsave_vint=181115.h5"
modal_params = map(x -> x.value, m.parameters) # The default parameters will be our "modal" parameters for this exercise

forecast_string = ""

# Set up regime switching and forecast settings
m <= Setting(:n_regimes, 2, true, "reg", "")   # How many regimes?
m <= Setting(:regime_switching, true)          # Need to set that the model does have regime-switching
m <= Setting(:forecast_block_size, 5000)
m <= Setting(:use_parallel_workers, true)
baseline_regime_dates = Dict{Int, Date}(1 => date_presample_start(m), 2 => Date(2020, 9, 30))
m <= Setting(:regime_dates, baseline_regime_dates)
setup_regime_switching_inds!(m)

# Some settings for alternative policy
if zlb_altpolicy
    # parametrize the alternative policy (AIT, in this case)
    m <= Setting(:pgap_value, 0.)
    m <= Setting(:pgap_type, :smooth_ait_gdp_alt)
    m <= Setting(:ygap_value, 12.)
    m <= Setting(:ygap_type, :smooth_ait_gdp_alt)
    m <= Setting(:smooth_ait_gdp_alt_ρ_smooth, 0.)
    m <= Setting(:ait_Thalf, 10.)
    m <= Setting(:gdp_Thalf, 10.)
    m <= Setting(:smooth_ait_gdp_alt_φ_y, 6.)
    m <= Setting(:smooth_ait_gdp_alt_φ_π, 6.)
end


for i in 1:2
    adj = (i == 2) ? .9 : 1. # adjustment to the size of the standard deviation

    # For each parameter, we need to instruct it that a regime exists.
    # Note that the regimes MUST be entered in order of the regimes
    # because regime values are stored in an OrderedDict, which is sorted
    # by insertion order.
    ModelConstructors.set_regime_val!(m[:σ_g], i, adj * m10[:σ_g].value)
    ModelConstructors.set_regime_val!(m[:σ_b], i, adj * m10[:σ_b].value)
    ModelConstructors.set_regime_val!(m[:σ_μ], i, adj * m10[:σ_μ].value)
    ModelConstructors.set_regime_val!(m[:σ_ztil], i, adj * m10[:σ_ztil].value)
    ModelConstructors.set_regime_val!(m[:σ_λ_f], i, adj * m10[:σ_λ_f].value)
    ModelConstructors.set_regime_val!(m[:σ_λ_w], i, adj * m10[:σ_λ_w].value)
    ModelConstructors.set_regime_val!(m[:σ_r_m], i, adj * m10[:σ_r_m].value)
    ModelConstructors.set_regime_val!(m[:σ_σ_ω], i, adj * m10[:σ_σ_ω].value)
    ModelConstructors.set_regime_val!(m[:σ_μ_e], i, adj * m10[:σ_μ_e].value)
    ModelConstructors.set_regime_val!(m[:σ_γ], i, adj * m10[:σ_γ].value)
    ModelConstructors.set_regime_val!(m[:σ_π_star], i, adj * m10[:σ_π_star].value)
    ModelConstructors.set_regime_val!(m[:σ_lr], i, adj * m10[:σ_lr].value)
    ModelConstructors.set_regime_val!(m[:σ_z_p], i, adj * m10[:σ_z_p].value)
    ModelConstructors.set_regime_val!(m[:σ_tfp], i, adj * m10[:σ_tfp].value)
    ModelConstructors.set_regime_val!(m[:σ_gdpdef], i, adj * m10[:σ_gdpdef].value)
    ModelConstructors.set_regime_val!(m[:σ_corepce], i, adj * m10[:σ_corepce].value)
    ModelConstructors.set_regime_val!(m[:σ_gdp], i, adj * m10[:σ_gdp].value)
    ModelConstructors.set_regime_val!(m[:σ_gdi], i, adj * m10[:σ_gdi].value)

    for j = 1:DSGE.n_mon_anticipated_shocks(m)
        ModelConstructors.set_regime_val!(m[Symbol("σ_r_m$(j)")], i, adj * m10[Symbol("σ_r_m$(j)")])
    end
end

df = load_data(m; check_empty_columns = false)

output_vars = [:histpseudo, :histobs, :histstdshocks,
               :hist4qpseudo, :hist4qobs, :histutpseudo,
               :forecastpseudo, :forecastobs, :forecastutpseudo,
               :forecast4qpseudo, :forecast4qobs, :forecaststdshocks]

if zlb_altpolicy
    Hbar = 4 # number of zlb regimes
    n_tempZLB_regimes = 3 + Hbar - 1 # 2 historical regimes, 4 quarters of zlb, and a switch to uncertain flexible ait in the last regime
    # Add parameter values for additional regimes
    for i in 3:(n_tempZLB_regimes+1)
        set_regime_val!(m[:σ_g], i, m10[:σ_g].value)
        set_regime_val!(m[:σ_b], i, m10[:σ_b].value)
        set_regime_val!(m[:σ_μ], i, m10[:σ_μ].value)
        set_regime_val!(m[:σ_ztil], i, m10[:σ_ztil].value)
        set_regime_val!(m[:σ_λ_f], i, m10[:σ_λ_f].value)
        set_regime_val!(m[:σ_λ_w], i, m10[:σ_λ_w].value)
        set_regime_val!(m[:σ_r_m], i, m10[:σ_r_m].value)
        set_regime_val!(m[:σ_σ_ω], i, m10[:σ_σ_ω].value)
        set_regime_val!(m[:σ_μ_e], i, m10[:σ_μ_e].value)
        set_regime_val!(m[:σ_γ], i, m10[:σ_γ].value)
        set_regime_val!(m[:σ_π_star], i, m10[:σ_π_star].value)
        set_regime_val!(m[:σ_lr], i, m10[:σ_lr].value)
        set_regime_val!(m[:σ_z_p], i, m10[:σ_z_p].value)
        set_regime_val!(m[:σ_tfp], i, m10[:σ_tfp].value)
        set_regime_val!(m[:σ_gdpdef], i, m10[:σ_gdpdef].value)
        set_regime_val!(m[:σ_corepce], i, m10[:σ_corepce].value)
        set_regime_val!(m[:σ_gdp], i, m10[:σ_gdp].value)
        set_regime_val!(m[:σ_gdi], i, m10[:σ_gdi].value)

        for j = 1:DSGE.n_mon_anticipated_shocks(m)
            set_regime_val!(m[Symbol("σ_r_m$(j)")], i, m10[Symbol("σ_r_m$(j)")].value)
        end
    end

    # then, set up new set of regime dates to implement a temporary zlb and uncertain policy change
    # Set up regime dates
    temp_regime_dates = deepcopy(get_setting(m, :regime_dates))
    temp_regime_dates[1] = date_presample_start(m)
    end_date = Date(2020, 12, 31) + Dates.Year(1) + Dates.Month(3)
    for (i, date) in zip(3:(n_tempZLB_regimes+1),
                         Date(2020, 12, 31):Dates.Month(3):end_date) # add forecast dates
        temp_regime_dates[i] = date
    end
    m <= Setting(:regime_dates, temp_regime_dates)
    m <= Setting(:regime_switching, true)
    setup_regime_switching_inds!(m)

    # Following two settings turn imperfect awareness
    m <= Setting(:uncertain_altpolicy, true) # turns imperfect awareness on
    m <= Setting(:uncertain_temporary_altpolicy, true) # needed to allow a temporary alternative policy with imperfect awareness

    # set the "alternative" rule for imperfect awareness
    # taylor_rule corresponds to the default monetary policy rule we use,
    # which has coefficients on the output gap, inflation gap, and output growth
    m <= Setting(:alternative_policies, AltPolicy[taylor_rule()])

    # Now set up settings for temp alt policy
    m <= Setting(:gensys2, true) # Temporary alternative policies use a special gensys algorithm
    m <= Setting(:replace_eqcond, true) # The gensys2 algo replaces eqcond matrices, so this step is required
    m <= Setting(:temporary_altpolicy_length, n_tempZLB_regimes - 2) # specifies number of ZLB regimes; required if there is
                                                               # further regime-switching after the ZLB ends (aside from
                                                               # the extra regime for the "lift-off" from ZLB)
    credvec = range(0., stop = 1., length = n_tempZLB_regimes - 2) # Credibility of ZLB increases as time goes on
    replace_eqcond = Dict{Int, EqcondEntry}() # Which eqcond to use in which periods
    for (i, reg) in enumerate(3:n_tempZLB_regimes)
        replace_eqcond[reg] = EqcondEntry(zero_rate(), [credvec[i], 1. - credvec[i]]) # Temp ZLB rule in these regimes
    end
    replace_eqcond[n_tempZLB_regimes + 1] = EqcondEntry(DSGE.smooth_ait_gdp_alt(), [1., 0.]) # switch to AIT in the final regime
    m <= Setting(:regime_eqcond_info, replace_eqcond) # Add mapping of regimes to new eqcond matrices

    # set up information sets
    m <= Setting(:tvis_information_set, vcat([1:1, 2:2],
                                             [i:get_setting(m, :n_regimes) for i in
                                              3:get_setting(m, :n_regimes)]))

    # pgap and ygap initialization
    df[!, :obs_pgap] .= NaN
    df[!, :obs_ygap] .= NaN
    ind_init = findfirst(df[!, :date] .== Date(2020, 6, 30))
    df[ind_init, :obs_pgap] = 0.125
    df[ind_init, :obs_ygap] = 12.

    # adjust nominal rates data during temporary ZLB
    tempzlb_quarters = DSGE.quarter_range(Date(2020, 12, 31),end_date)
    start_ind = findfirst(df[!, :date] .== Date(2020, 12, 31))
    if !isnothing(start_ind)
        inds_tempzlb = start_ind:findfirst(date_forecast_start(m))
        df[inds_tempzlb, :obs_nominalrate] .= NaN
        df[inds_tempzlb, [Symbol("obs_nominalrate$i") for i in 1:n_mon_anticipated_shocks(m)]] .= NaN
    end

    fcast_tempzlb = DSGE.forecast_one_draw(m, :mode, :none, output_vars, modal_params,
                                           df; regime_switching = true, n_regimes = get_setting(m, :n_regimes))
end

if make_plots
    plots_dict = Dict()

    # Define some plotting functions
    function make_comp_plot(m::AbstractDSGEModel, plot_mat::OrderedDict{String, Dict{Symbol, Array{Float64}}},
                            obs::Symbol, product::Symbol,
                            dates_plot, data_plot,
                            fcast_date, nfcast;
                            styles::Vector{Symbol} = repeat([:solid], length(keys(plot_mat))),
                            colors::Vector{Symbol} = repeat([:red], length(keys(plot_mat))))
        series_ind = if product ==:forecastobs
            m.observables[obs]
        else
            m.pseudo_observables[obs]
        end

        if product == :forecastobs
            dates_plot = data_plot[!, :date]
            histobs_plot = obs == :obs_hours ? data_plot[!, obs] : 4 .* data_plot[!, obs]
            fcast_date_tmp = vcat(dates_plot[end - 1], fcast_date)
            dates_plot = map(x -> DSGE.prev_quarter(x), dates_plot)
            p = Plots.plot(dates_plot, histobs_plot, color = :black,
                           label = "Data", linewidth = 3, legend = :bottomright, left_margin = 20px)
        else
            fcast_date_tmp = fcast_date
            p = Plots.plot()
        end

        for (i, (key, val)) in enumerate(plot_mat)

            if obs == :obs_hours || product == :forecastpseudo
                forecast_plot = val[product][series_ind, :][1:nfcast]
            else
                forecast_plot = 4 .* val[product][series_ind, :][1:nfcast]
            end
            if product == :forecastobs
                forecast_plot = vcat(histobs_plot[end], forecast_plot)
            end

            plot!(fcast_date_tmp, forecast_plot, color = colors[i],
                  label = key, linewidth = 1, linestyle = styles[i])
            if obs == :obs_nominalrate
                ylims!((-2, 4))
            end
        end
        return p
    end

    data_plot  = df[df[:, :date] .>= Date(2016, 3, 31), :]
    fcast_date = map(x -> DSGE.prev_quarter(x), DSGE.get_quarter_ends(date_forecast_start(m), DSGE.quartertodate("2031-Q4")))
    nfcast     = length(fcast_date)
    dates_plot = map(x -> DSGE.prev_quarter(x), data_plot[!, :date])

    for obs in [:obs_gdp, :obs_corepce, :obs_gdpdeflator, :obs_nominalrate, :obs_investment, :obs_consumption,
                :obs_spread, :obs_wages, :obs_hours]
        plots_dict[obs] = make_comp_plot(m, OrderedDict("$(pol_str)" => fcast_tempzlb),
                                         obs, :forecastobs,
                                         dates_plot, data_plot, fcast_date, nfcast,
                                         colors = [:red, :blue, :purple])
    end
end

if save_plots
    if !isdir(figurespath(m, "forecast"))
        mkdir(figurespath(m, "forecast"))
    end

    for (k, p) in plots_dict
        Plots.savefig(p, joinpath(figuespath(m, "forecast"), "$(pol_str)_$(string(k)).pdf"))
    end
end

nothing
