using DSGE, ModelConstructors, Dates, CSV, DataFrames, Plots, Test
include("tvcred_parameterize.jl")

# Forecast settings
fcast_date = DSGE.quartertodate("2020-Q4")
date_fcast_end  = iterate_quarters(fcast_date, 60)

# Flexible AIT Rule
imperfect_cred_flexait = .33
imperfect_cred_taylor  = 1. - imperfect_cred_flexait
φ_π                    = 10.
φ_y                    = 6.
Thalf                  = 10.
ρ_smooth               = 0.7
pgap_init              = 0.125
ygap_init              = 12.
pgap_ygap_init_date    = Date(2020, 6, 30)
start_zlb_date         = Date(2021, 3, 31)
end_zlb_date           = Date(2022, 12, 31)

# additional settings to implement Flexible AIT rule
flexait_custom = Dict{Symbol, Setting}(:add_initialize_pgap_ygap_pseudoobs => Setting(:add_initialize_pgap_ygap_pseudoobs, true),
                                       :add_pgap => Setting(:add_pgap, true), :add_ygap => Setting(:add_ygap, true),
                                       :add_altpolicy_pgap => Setting(:add_altpolicy_pgap, true),
                                       :add_altpolicy_ygap => Setting(:add_altpolicy_ygap, true),
                                       :forecast_horizons => Setting(:forecast_horizons,
                                                                     subtract_quarters(date_fcast_end, fcast_date)),
                                       :forecast_smoother => Setting(:forecast_smoother, :carter_kohn),
                                       :contemporaneous_and_proportional_antshocks =>
                                       Setting(:contemporaneous_and_proportional_antshocks, Symbol[:biidc]),
                                       :antshocks =>
                                       Setting(:antshocks, Dict{Symbol, Int}(:biidc => 1,
                                                                             :φ => 1, :ziid => 1)),
                                       :ant_eq_mapping =>
                                       Setting(:ant_eq_mapping, Dict{Symbol, Symbol}(:biidc => :biidc,
                                                                                     :φ => :φ,
                                                                                     :ziid => :ziid)),
                                       :ant_eq_E_mapping =>
                                       Setting(:ant_eq_E_mapping, Dict{Symbol, Symbol}(:φ => :Eφ)),
                                       :proportional_antshocks =>
                                       Setting(:proportional_antshocks, [:biidc, :φ, :ziid]),
                                       :n_anticipated_obs_gdp => Setting(:n_anticipated_obs_gdp, 1),
                                       :add_anticipated_obs_gdp => Setting(:add_anticipated_obs_gdp, true),
                                       :flexible_ait_2020Q3_policy_change =>
                                       Setting(:flexible_ait_2020Q3_policy_change, false))

# Initialize model object
m = Model1002("ss59"; custom_settings = flexait_custom)
usual_model_settings!(m, "200001", cdvt = "200001", fcast_date = fcast_date)
m <= Setting(:time_varying_trends, true)
get_setting(m, :regime_dates)[5] = Date(2020, 12, 31)
setup_regime_switching_inds!(m, cond_type = :full)

# Parameterize
tvcred_parameterize!(m)

# Load data
df = CSV.read(joinpath(dirname(@__FILE__), "..", "reference", "timevaryingcred_df.csv"), DataFrame)
df[!, :obs_pgap] .= NaN # NaN out just to be safe
df[!, :obs_ygap] .= NaN
ind_init = findfirst(df[!, :date] .== pgap_ygap_init_date)

# Set up Flexible AIT with a temporary ZLB. Note that df will be altered
θ = (φ_π = φ_π, φ_y = φ_y, Thalf = Thalf, pgap = pgap_init,
     ygap = ygap_init, ρ_smooth = ρ_smooth,
     cred = imperfect_cred_flexait,
     historical_policy = DSGE.taylor_rule()) # NamedTuple with parameters of the Flexible AIT policy
df[ind_init, :obs_pgap] = -θ[:pgap] # initialize pgap and ygap in the date
df[ind_init, :obs_ygap] = -θ[:ygap] # pgap_ygap_init_date, default to 2020:Q2

m <= Setting(:flexible_ait_φ_π, θ[:φ_π])
m <= Setting(:flexible_ait_φ_π, θ[:φ_y])
m <= Setting(:ait_Thalf, θ[:Thalf])
m <= Setting(:gdp_Thalf, θ[:Thalf])
m <= Setting(:pgap_value, θ[:pgap])
m <= Setting(:pgap_type, :flexible_ait)
m <= Setting(:ygap_value, θ[:ygap])
m <= Setting(:ygap_type, :flexible_ait)
m <= Setting(:flexible_ait_ρ_smooth, θ[:ρ_smooth])
m <= Setting(:alternative_policy, DSGE.flexible_ait())
m <= Setting(:alternative_policy_weights, [θ[:cred], 1. - θ[:cred]])
m <= Setting(:alternative_policies, AltPolicy[θ[:historical_policy]])
m <= Setting(:skip_altpolicy_state_init, true)

## Set up temporary ZLB
m <= Setting(:gensys2, true)

# Set up credibility for ZLB and alternative policy, if applicable
m <= Setting(:uncertain_zlb, true)
m <= Setting(:uncertain_altpolicy, true)
gensys2_first_regime = get_setting(m, :n_regimes) + 1
get_setting(m, :regime_dates)[gensys2_first_regime] = start_zlb_date
n_zlb_reg = DSGE.subtract_quarters(end_zlb_date, start_zlb_date) + 1

# Set up replace_eqcond_func_dict
m <= Setting(:replace_eqcond, true)
reg_dates = deepcopy(get_setting(m, :regime_dates))
replace_eqcond_func_dict = Dict{Int, Function}()
for (regind, date) in zip(gensys2_first_regime:(n_zlb_reg - 1 + gensys2_first_regime), # See comments starting at line 57
                          DSGE.quarter_range(reg_dates[gensys2_first_regime],
                                             DSGE.iterate_quarters(reg_dates[gensys2_first_regime], n_zlb_reg - 1)))
    reg_dates[regind] = date
    replace_eqcond_func_dict[regind] = zero_rate_replace_eq_entries
end
reg_dates[n_zlb_reg + gensys2_first_regime] = DSGE.iterate_quarters(reg_dates[gensys2_first_regime], n_zlb_reg)
replace_eqcond_func_dict[n_zlb_reg + gensys2_first_regime] = DSGE.flexible_ait_replace_eq_entries
nreg0 = length(reg_dates)
m <= Setting(:regime_dates,             reg_dates)
m <= Setting(:replace_eqcond_func_dict, replace_eqcond_func_dict)

# Set up regime indices and extra regimes for parameters
setup_regime_switching_inds!(m; cond_type = :full)
set_regime_vals_fnct(m, get_setting(m, :n_regimes))

# Set up TVIS information set
m <= Setting(:tvis_information_set, vcat([i:i for i in 1:(gensys2_first_regime - 1)],
                                         [i:get_setting(m, :n_regimes) for i in
                                          gensys2_first_regime:get_setting(m, :n_regimes)]))
# CHANGING THE LOCATION OF THE LINE SETTING TEMPORARY ZLB LENGTH CHANGES THE FORECAST SO THERE'S A BUG SOMEWHERE
# IF TEMPORARY ZLB LENGTH ISN'T SET FOR THE FIXED CRED FORECAST, YOU GET DIFFERENT RESULTS
# m <= Setting(:temporary_zlb_length, n_zlb_reg)

output_vars = [:histobs, :histpseudo, :forecastobs, :forecastpseudo]
modal_params = map(x -> x.value, m.parameters)
m <= Setting(:alternative_policy_weights, [0., 1.])
@assert !haskey(m.settings, :alternative_policy_varying_weights)
outp0 = out = DSGE.forecast_one_draw(m, :mode, :full, output_vars, modal_params, df,
                             regime_switching = true, n_regimes = get_setting(m, :n_regimes))
m <= Setting(:alternative_policy_weights, [θ[:cred], 1. - θ[:cred]])
outp33 = out = DSGE.forecast_one_draw(m, :mode, :full, output_vars, modal_params, df,
                             regime_switching = true, n_regimes = get_setting(m, :n_regimes))

for i in (nreg0 - 1):(nreg0 + 15)
    reg_dates[i] = DSGE.iterate_quarters(reg_dates[nreg0], i - nreg0)
    replace_eqcond_func_dict[i] = DSGE.flexible_ait_replace_eq_entries
end
m <= Setting(:regime_dates,             reg_dates)
m <= Setting(:replace_eqcond_func_dict, replace_eqcond_func_dict)

m <= Setting(:temporary_zlb_length, n_zlb_reg) # CHANGING THE LOCATION OF TEMPORARY ZLB LENGTH CHANGES THE FORECAST
m <= Setting(:alternative_policy_varying_weights,
             Dict(k => [.33, 1. - .33] for k in keys(get_setting(m, :replace_eqcond_func_dict))))
outp33_tv = DSGE.forecast_one_draw(m, :mode, :full, output_vars, modal_params, df,
                                 regime_switching = true, n_regimes = get_setting(m, :n_regimes))
for k in keys(outp33)
#    @test outp33[k] ≈ outp33_tv[k]
end

m <= Setting(:alternative_policy_varying_weights,
             Dict(k => [0., 1.] for k in keys(get_setting(m, :replace_eqcond_func_dict))))
credvec = collect(range(0., stop = 1., length = 17))
for (i, k) in enumerate(sort!(collect(keys(replace_eqcond_func_dict))))
    if (get_setting(m, :temporary_zlb_length) - 1) < i
        get_setting(m, :alternative_policy_varying_weights)[k] = [credvec[i - (get_setting(m, :temporary_zlb_length) - 1)], 1. - credvec[i - (get_setting(m, :temporary_zlb_length) - 1)]]
    end
end
out1 = DSGE.forecast_one_draw(m, :mode, :full, output_vars, modal_params, df,
                              regime_switching = true, n_regimes = get_setting(m, :n_regimes))

m <= Setting(:alternative_policy_varying_weights,
             Dict(k => [0., 1.] for k in keys(get_setting(m, :replace_eqcond_func_dict))))
credvec = collect(range(0., stop = .5, length = 17))
for (i, k) in enumerate(sort!(collect(keys(replace_eqcond_func_dict))))
    if (get_setting(m, :temporary_zlb_length) - 1) < i
        get_setting(m, :alternative_policy_varying_weights)[k] = [credvec[i - (get_setting(m, :temporary_zlb_length) - 1)], 1. - credvec[i - (get_setting(m, :temporary_zlb_length) - 1)]]
    end
end
out2 = DSGE.forecast_one_draw(m, :mode, :full, output_vars, modal_params, df,
                              regime_switching = true, n_regimes = get_setting(m, :n_regimes))

m <= Setting(:alternative_policy_varying_weights,
             Dict(k => [0., 1.] for k in keys(get_setting(m, :replace_eqcond_func_dict))))
credvec = collect(range(0., stop = .25, length = 17))
for (i, k) in enumerate(sort!(collect(keys(replace_eqcond_func_dict))))
    if (get_setting(m, :temporary_zlb_length) - 1) < i
        get_setting(m, :alternative_policy_varying_weights)[k] = [credvec[i - (get_setting(m, :temporary_zlb_length) - 1)], 1. - credvec[i - (get_setting(m, :temporary_zlb_length) - 1)]]
    end
end
out3 = DSGE.forecast_one_draw(m, :mode, :full, output_vars, modal_params, df,
                              regime_switching = true, n_regimes = get_setting(m, :n_regimes))

plot_dict = Dict()
for k in [:obs_gdp, :obs_corepce, :obs_nominalrate, :obs_pgap, :obs_ygap]
    plot_dict[k] = plot()
    plot!(1:30, outp0[:forecastobs][m.observables[k], 1:30], label = "Fixed Cred = 0%", linewidth = 2,
          legend = :topright)
    plot!(1:30, outp33[:forecastobs][m.observables[k], 1:30], label = "System", linewidth = 2)
    plot!(1:30, out1[:forecastobs][m.observables[k], 1:30], label = "TV Cred to 100%", linewidth = 2)
    plot!(1:30, out2[:forecastobs][m.observables[k], 1:30], label = "TV Cred to 50%", linewidth = 2)
    plot!(1:30, out3[:forecastobs][m.observables[k], 1:30], label = "TV Cred to 25%", linewidth = 2)
end
for k in [:NaturalRate, :OutputGap]
    plot_dict[k] = plot()
    plot!(1:30, outp0[:forecastpseudo][m.pseudo_observables[k], 1:30], label = "Fixed Cred = 0%", linewidth = 2,
          legend = :bottomright)
    plot!(1:30, outp33[:forecastpseudo][m.pseudo_observables[k], 1:30], label = "System", linewidth = 2)
    plot!(1:30, out1[:forecastpseudo][m.pseudo_observables[k], 1:30], label = "TV Cred to 100%", linewidth = 2)
    plot!(1:30, out2[:forecastpseudo][m.pseudo_observables[k], 1:30], label = "TV Cred to 50%", linewidth = 2)
    plot!(1:30, out3[:forecastpseudo][m.pseudo_observables[k], 1:30], label = "TV Cred to 25%", linewidth = 2)
end

#=
mkdir("tvcred_figs")
for (k, v) in plot_dict
    savefig(v, "tvcred_figs/$(k).pdf")
end
=#
nothing
