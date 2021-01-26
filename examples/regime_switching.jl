using ClusterManagers, DSGE, ModelConstructors, FileIO, Plots, Dates, OrderedCollections, Distributed
using Plots.PlotMeasures
GR.inline("pdf")
fn = dirname(@__FILE__)

# This script shows how to use exogenous regime-switching
# for forecasting as well as permanent and
# temporary alternative policies. The format
# of this script mirrors make_packet.jl,
# so users should familiarize themselves with that
# script first.

# In this script, we consider a case where
# in 1990:Q1, the standard deviations of
# all shocks are set to 0.9 of their values in
# the period 1959:Q3-1989:Q4, but by 2015:Q1,
# the standard deviations return to the pre-1990 values.
# Furthermore, we consider the effects of a permanent and temporary
# alternative policy of nominal GDP (NGDP) targeting.

# What do you want to do?
run_full_forecast = true  # Run full distribution forecast
temp_alt          = true  # Temporary NGDP targeting?
perm_alt          = true  # Permanent NGDP targeting?
make_plots        = true  # Make plots
save_plots        = false # Save plots to figurespath(m, "forecast")
add_workers       = false # Run in parallel
n_workers         = 10

# Initialize model objects and desired settings
custom_settings = Dict{Symbol, Setting}(:add_altpolicy_pgap => Setting(:add_altpolicy_pgap, true))
m = Model1002("ss10"; custom_settings = custom_settings)   # We will directly construct the matrix of parameter draws
m10 = Model1002("ss10") # b/c loading from a saved estimation file has not been fully implemented

for model in [m, m10]
    # Not strictly necessary to have the same settings, but it makes it easier for us
    model <= Setting(:sampling_method, :MH)
    usual_model_settings!(model, "181115"; cdvt = "181115", fcast_date = quartertodate("2018-Q4"))
    model <= Setting(:use_population_forecast, false)
    model <= Setting(:saveroot, joinpath("$(fn)", "..", "save"))
    model <= Setting(:dataroot, joinpath("$(fn)", "..", "save", "input_data"))
    model <= Setting(:date_forecast_end, quartertodate("2030-Q1"))
    model <= Setting(:forecast_horizons, DSGE.subtract_quarters(date_forecast_end(model), date_forecast_start(model)) + 1)
    model <= Setting(:forecast_block_size, 40)
    if get_setting(model, :sampling_method) == :SMC
        model <= Setting(:forecast_jstep, 1)
    end
end

forecast_string = ""

# Set up regime switching and forecast settings
m <= Setting(:n_regimes, 3, true, "reg", "")   # How many regimes?
m <= Setting(:regime_switching, true)          # Need to set that the model does have regime-switching
# m <= Setting(:regime_switching_ndraws, 100000) # b/c we don't have an estimation file, need this setting for dimension purposes
m <= Setting(:forecast_block_size, 5000)
m <= Setting(:use_parallel_workers, true)
baseline_regime_dates = Dict{Int, Date}(1 => date_presample_start(m), 2 => Date(1990, 3, 31),
                                        3 => Date(2015, 3, 31))
m <= Setting(:regime_dates, baseline_regime_dates)
setup_regime_switching_inds!(m)

# Some settings for alternative policy
if temp_alt || perm_alt
    # Define values for creating the alternative policy, for either case
    policy               = :ngdp
    replace_policy       = DSGE.ngdp_replace_eq_entries
    policy_eqcond        = DSGE.ngdp_eqcond
    policy_solve         = DSGE.ngdp_solve
    policy_forecast_init = DSGE.ngdp_forecast_init
    pol_str              = "NGDP"
end

# Set up parameters for regime-switching.
# See the full forecast code block for the construction
# of the matrix of parameter draws
overrides = forecast_input_file_overrides(m10)
overrides[:full] = "$(fn)/../test/reference/mhsave_vint=181115.h5"
modal_params = map(x -> x.value, m10.parameters) # The default parameters will be our "modal" parameters for this exercise
θ10 = load_draws(m10, :full)

para_regs = perm_alt ? (1:4) : (1:3) # if permanent altpolicy occurs, then we need to add a fourth regime
for i in para_regs
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

# Full-distribution forecast
if run_full_forecast

    # Set up matrix of parameter draws
    m <= Setting(:forecast_jstep, 1)
    m <= Setting(:forecast_block_size, 100)
    param_mat = zeros(100, ModelConstructors.n_parameters_regime_switching(m)) # Need to use a special method to count the additional regimes
    n_pvec = n_parameters(m) # This is just length(m.parameters), so it counts the "unique" number of parameters, excluding regimes
    for i in 1:size(param_mat, 1)
        # Start w/regime 1 b/c all regime 1 parameters are grouped together
        param_mat[i, 1:n_pvec] = map(x -> (abs(x) < 1e-15) ? 0. : x, θ10[i, :])

        # Now do the other regimes, which grouped by parameter rather than by regime
        j = n_pvec # To help index for regime-switching parameters
        for para in m.parameters
            if !isempty(para.regimes)
                for (k, v) in para.regimes[:value] # If there exists multiple regimes
                    if k > 1
                        j += 1
                        drawv = θ10[i, m10.keys[para.key]] # Find the column for the parameter in θ10 using the parameter's key
                        param_mat[i, j] = (k  == 2) ? .9 * drawv : drawv # Adjust the value
                        if abs(param_mat[i, j]) < 1e-15 # Make sure essentially zero parameters are actually zero
                            param_mat[i, j] = 0.
                        end
                    end
                end
            end
        end
    end

    output_vars = [:histpseudo, :histobs, :histstdshocks,
                   :hist4qpseudo, :hist4qobs, :histutpseudo,
                   :forecastpseudo, :forecastobs, :forecastutpseudo,
                   :forecast4qpseudo, :forecast4qobs, :forecaststdshocks]

    if add_workers
        my_procs = addprocs(n_workers)
        @everywhere using DSGE, OrderedCollections
    end

    usual_model_forecast(m, :full, :none, output_vars,
                         forecast_string = forecast_string,
                         density_bands = [.5, .6, .68, .7, .8, .9],
                         check_empty_columns = false,
                         params = param_mat) # Need to pass in parameters directly, see ?forecast_one

    if add_workers
        rmprocs(my_procs)
    end

    # For comparison to the alternative policies. We want the raw output here, so we call
    # forecast_one_draw rataher than the wrappers usual_model_forecast or forecast_one
    # Note m is already calibrated to the "mode", so passing in parameters of same length
    # as m.parameters will just trigger the non-regime-switching updating of the parameters in m.
    # If the parameters of m were not already calibrated, then modal_params would need to be extended
    # to include the correct regime-switching values (rather than just map(x -> x.value, m10.parameter))
    fcast = DSGE.forecast_one_draw(m, :mode, :none, output_vars, modal_params,
                                   df; regime_switching = true, n_regimes = get_setting(m, :n_regimes))
end

if perm_alt
    # m <= Setting(:replace_eqcond, false)  # Make sure this setting is false
    m <= Setting(:pgap_type, policy)      # Use NGDP
    m <= Setting(:pgap_value, 12.0)       # parametrizing the NGDP rule
    m <= Setting(:regime_switching, true) # Regime-switching should still be on b/c historical regimes
    m <= Setting(:regime_dates, baseline_regime_dates) # Make sure to use historical regimes
    setup_regime_switching_inds!(m)

    m <= Setting(:gensys2, false) # Make sure we use the normal gensys algorithm

    setup_permanent_altpol!(m, AltPolicy(policy, policy_eqcond, policy_solve,
                                         forecast_init = policy_forecast_init)) # Set up permanent alternative policy
    fcast_permalt = DSGE.forecast_one_draw(m, :mode, :none, output_vars, modal_params,
                                           df; regime_switching = true, n_regimes = get_setting(m, :n_regimes))
end

if temp_alt
    # First, need to set up new set of regime dates to implement a temporary alternative policy (temp alt policy)
    temp_n_regimes = 16 # 3 historical regimes, 12 quarters or 3 years of the temp alt policy, and a return to normal in the last eregime
    temp_regime_dates = Dict{Int, Date}()
    temp_regime_dates[1] = date_presample_start(m)
    for (k, v) in baseline_regime_dates # Add the historical regimes
        temp_regime_dates[k] = v
    end
    end_date = Date(2018, 12, 31) + Dates.Year(3) + Dates.Month(3)
    for (i, date) in zip(4:temp_n_regimes,
                         Date(2018, 12, 31):Dates.Month(3):end_date) # add forecast dates
        temp_regime_dates[i] = date
    end
    m <= Setting(:regime_dates, temp_regime_dates)
    m <= Setting(:regime_switching, true)
    setup_regime_switching_inds!(m)

    # Add parameter values for additional regimes
    for i in 4:temp_n_regimes
        ModelConstructors.set_regime_val!(m[:σ_g], i, m10[:σ_g].value)
        ModelConstructors.set_regime_val!(m[:σ_b], i, m10[:σ_b].value)
        ModelConstructors.set_regime_val!(m[:σ_μ], i, m10[:σ_μ].value)
        ModelConstructors.set_regime_val!(m[:σ_ztil], i, m10[:σ_ztil].value)
        ModelConstructors.set_regime_val!(m[:σ_λ_f], i, m10[:σ_λ_f].value)
        ModelConstructors.set_regime_val!(m[:σ_λ_w], i, m10[:σ_λ_w].value)
        ModelConstructors.set_regime_val!(m[:σ_r_m], i, m10[:σ_r_m].value)
        ModelConstructors.set_regime_val!(m[:σ_σ_ω], i, m10[:σ_σ_ω].value)
        ModelConstructors.set_regime_val!(m[:σ_μ_e], i, m10[:σ_μ_e].value)
        ModelConstructors.set_regime_val!(m[:σ_γ], i, m10[:σ_γ].value)
        ModelConstructors.set_regime_val!(m[:σ_π_star], i, m10[:σ_π_star].value)
        ModelConstructors.set_regime_val!(m[:σ_lr], i, m10[:σ_lr].value)
        ModelConstructors.set_regime_val!(m[:σ_z_p], i, m10[:σ_z_p].value)
        ModelConstructors.set_regime_val!(m[:σ_tfp], i, m10[:σ_tfp].value)
        ModelConstructors.set_regime_val!(m[:σ_gdpdef], i, m10[:σ_gdpdef].value)
        ModelConstructors.set_regime_val!(m[:σ_corepce], i, m10[:σ_corepce].value)
        ModelConstructors.set_regime_val!(m[:σ_gdp], i, m10[:σ_gdp].value)
        ModelConstructors.set_regime_val!(m[:σ_gdi], i, m10[:σ_gdi].value)

        for j = 1:DSGE.n_mon_anticipated_shocks(m)
            ModelConstructors.set_regime_val!(m[Symbol("σ_r_m$(j)")], i, m10[Symbol("σ_r_m$(j)")].value)
        end
    end

    # Now set up settings for temp alt policy
    m <= Setting(:gensys2, true) # Temporary alternative policies use a special gensys algorithm
    m <= Setting(:replace_eqcond, true) # This new gensys algo replaces eqcond matrices, so this step is required
    regime_eqcond_info = Dict{Int, DSGE.EqcondEntry}() # Which eqcond to use in which periods
    for i in 4:(temp_n_regimes - 1)
        regime_eqcond_info[i] = DSGE.EqcondEntry(DSGE.ngdp()) # Use the alternative policy for these regimes!
    end
    regime_eqcond_info[temp_n_regimes] = DSGE.EqcondEntry(default_policy()) # Use the alternative policy for these regimes!
    m <= Setting(:regime_eqcond_info, regime_eqcond_info) # Add mapping of regimes to new eqcond matrices

    m <= Setting(:pgap_value, 12.0)  # parametrizing the NGDP rule
    m <= Setting(:pgap_type, policy)

    fcast_tempalt = DSGE.forecast_one_draw(m, :mode, :none, output_vars, modal_params,
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

    function make_ngdp_plot(m::AbstractDSGEModel, plot_mat::OrderedDict{String, Dict{Symbol, Array{Float64}}},
                            dates_plot, data_plot,
                            fcast_date, nfcast;
                            styles::Vector{Symbol} = repeat([:solid], length(keys(plot_mat))),
                            colors::Vector{Symbol} = repeat([:red], length(keys(plot_mat))))
        fcast_permalt_date = map(x -> DSGE.prev_quarter(x), DSGE.get_quarter_ends(date_forecast_start(m),
                                                                                  DSGE.quartertodate("2024-Q4")))
        nfcast_permalt     = length(fcast_permalt_date) - 1

        p = Plots.plot()
        for (i, (key, val)) in enumerate(plot_mat)
            plotobsgdp = vcat(df[end - 1, :obs_gdp], val[:forecastobs][m.observables[:obs_gdp], 1:nfcast_permalt])
            plotpi = vcat(df[end - 1, :obs_gdpdeflator],
                          val[:forecastobs][m.observables[:obs_gdpdeflator], 1:nfcast_permalt])
            ngdp_level = cumprod(exp.(plotobsgdp ./ 100.) + plotpi ./ 100.)
            ngdp_level ./= ngdp_level[1]
            ngdp_level .-= 1.
            ngdp_level .*= 100

            plot!(fcast_permalt_date[1:end], ngdp_level[1:end], color = colors[i], label = key, linewidth = 3, linestyle = styles[i])
        end
        return p
    end

    data_plot  = df[df[:, :date] .>= Date(2016, 3, 31), :]
    fcast_date = map(x -> DSGE.prev_quarter(x), DSGE.get_quarter_ends(date_forecast_start(m), date_forecast_end(m)))
    nfcast     = length(fcast_date)
    dates_plot = map(x -> DSGE.prev_quarter(x), data_plot[!, :date])

    for obs in [:obs_gdp, :obs_corepce, :obs_gdpdeflator, :obs_nominalrate, :obs_investment, :obs_consumption,
                :obs_spread, :obs_wages, :obs_hours]
        plots_dict[obs] = make_comp_plot(m, OrderedDict("Temp $(pol_str)" => fcast_tempalt,
                                                          "Baseline" => fcast,
                                                          "$(pol_str)" => fcast_permalt),
                                         obs, :forecastobs,
                                         dates_plot, data_plot, fcast_date, nfcast,
                                         colors = [:red, :blue, :purple])
    end

    plots_dict[:cum_NominalGDPLevel] = make_ngdp_plot(m, OrderedDict("Temp $(pol_str)" => fcast_tempalt,
                                                                       "Baseline" => fcast,
                                                                       "$(pol_str)" => fcast_permalt),
                                                      dates_plot, data_plot, fcast_date, nfcast,
                                                      colors = [:red, :blue, :purple])
end

if save_plots
    if !isdir(figurespath(m, "forecast"))
        mkdir(figurespath(m, "forecast"))
    end

    for (k, p) in plots_dict
        Plots.savefig(p, joinpath(figurespath(m, "forecast"), "$(pol_str)_$(string(k)).pdf"))
    end
end


nothing
