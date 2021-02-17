using DSGE, ModelConstructors, Dates, CSV, DataFrames, Plots, Test, HDF5
include("tvcred_parameterize.jl")

regenerate_reference_forecasts = false

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
                                       :forecast_smoother => Setting(:forecast_smoother, :durbin_koopman),
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
                                        :meas_err_anticipated_obs_gdp =>
                                        Setting(:meas_err_anticipated_obs_gdp, 1.),
                                       :flexible_ait_policy_change =>
                                       Setting(:flexible_ait_policy_change, false))

# Initialize model object
m = Model1002("ss59"; custom_settings = flexait_custom)
usual_model_settings!(m, "200001", cdvt = "200001", fcast_date = fcast_date)
m <= Setting(:time_varying_trends, true)
get_setting(m, :regime_dates)[5] = Date(2020, 12, 31)
setup_regime_switching_inds!(m, cond_type = :full)

m10 = Model1002("ss10") # for help initializing parameters of m
θ10 = h5read(joinpath(dirname(@__FILE__), "..", "reference", "tvcred_reference_forecast.h5"), "para")
DSGE.update!(m10, θ10)
update_vals = zeros(Float64, length(m.parameters))
for k in map(x-> x.key, m10.parameters)
    update_vals[m.keys[k]] = m10[k].value
end
DSGE.update!(m, update_vals) # make sure parameters match Model 1002 ss10, w/zeros for other parameters

# Parameterize
tvcred_parameterize!(m)

# Load data
df = CSV.read(joinpath(dirname(@__FILE__), "..", "reference", "timevaryingcred_df.csv"), DataFrame)
df[!, :obs_pgap] .= NaN # NaN out just to be safe
df[!, :obs_ygap] .= NaN
ind_init = findfirst(df[!, :date] .== pgap_ygap_init_date)
df[end, :obs_gdp] = df[end - 1, :obs_gdp1] # To match reference forecast

# Set up Flexible AIT with a temporary ZLB. Note that df will be altered
θ = (φ_π = φ_π, φ_y = φ_y, Thalf = Thalf, pgap = pgap_init,
     ygap = ygap_init, ρ_smooth = ρ_smooth,
     cred = imperfect_cred_flexait,
     historical_policy = DSGE.taylor_rule()) # NamedTuple with parameters of the Flexible AIT policy
df[ind_init, :obs_pgap] = -θ[:pgap] # initialize pgap and ygap in the date
df[ind_init, :obs_ygap] = -θ[:ygap] # pgap_ygap_init_date, default to 2020:Q2

m <= Setting(:flexible_ait_φ_π, θ[:φ_π])
m <= Setting(:flexible_ait_φ_y, θ[:φ_y])
m <= Setting(:ait_Thalf, θ[:Thalf])
m <= Setting(:gdp_Thalf, θ[:Thalf])
m <= Setting(:pgap_value, θ[:pgap])
m <= Setting(:pgap_type, :flexible_ait)
m <= Setting(:ygap_value, θ[:ygap])
m <= Setting(:ygap_type, :flexible_ait)
m <= Setting(:flexible_ait_ρ_smooth, θ[:ρ_smooth])
m <= Setting(:alternative_policies, AltPolicy[θ[:historical_policy]])
m <= Setting(:skip_altpolicy_state_init, true)

## Set up temporary ZLB
m <= Setting(:gensys2, true)
m <= Setting(:uncertain_temp_altpol, true)
m <= Setting(:uncertain_altpolicy, true)

# Set regime_dates and regime_eqcond_info for temp ZLB to flexible AIT
gensys2_first_regime = get_setting(m, :n_regimes) + 1
get_setting(m, :regime_dates)[gensys2_first_regime] = start_zlb_date
n_zlb_reg = DSGE.subtract_quarters(end_zlb_date, start_zlb_date) + 1
m <= Setting(:replace_eqcond, true)
reg_dates = deepcopy(get_setting(m, :regime_dates))
regime_eqcond_info = Dict{Int, DSGE.EqcondEntry}()
for (regind, date) in zip(gensys2_first_regime:(n_zlb_reg - 1 + gensys2_first_regime), # See comments starting at line 57
                          DSGE.quarter_range(reg_dates[gensys2_first_regime],
                                             DSGE.iterate_quarters(reg_dates[gensys2_first_regime], n_zlb_reg - 1)))
    reg_dates[regind] = date
    regime_eqcond_info[regind] = DSGE.EqcondEntry(DSGE.zero_rate(), [0.5, 0.5])
end
reg_dates[n_zlb_reg + gensys2_first_regime] = DSGE.iterate_quarters(reg_dates[gensys2_first_regime], n_zlb_reg)
regime_eqcond_info[n_zlb_reg + gensys2_first_regime] = DSGE.EqcondEntry(DSGE.flexible_ait(), [0.5, 0.5])
nreg0 = length(reg_dates)
m <= Setting(:regime_dates,             reg_dates)
m <= Setting(:regime_eqcond_info, regime_eqcond_info)
setup_regime_switching_inds!(m; cond_type = :full)
set_regime_vals_fnct(m, get_setting(m, :n_regimes))
m <= Setting(:temporary_altpol_length, n_zlb_reg)

# Now add additional regimes of flexible AIT to allow time-varying credibility
for i in nreg0:(nreg0 + 15)
    reg_dates[i] = DSGE.iterate_quarters(reg_dates[nreg0], i - nreg0)
    regime_eqcond_info[i] = DSGE.EqcondEntry(DSGE.flexible_ait(), [0.5, 0.5])
end
m <= Setting(:regime_dates,             reg_dates)
m <= Setting(:regime_eqcond_info, regime_eqcond_info)
setup_regime_switching_inds!(m; cond_type = :full)
set_regime_vals_fnct(m, get_setting(m, :n_regimes))

# Set up TVIS information set
m <= Setting(:tvis_information_set, vcat([i:i for i in 1:(gensys2_first_regime - 1)],
                                         [i:get_setting(m, :n_regimes) for i in
                                          gensys2_first_regime:get_setting(m, :n_regimes)]))

# Run forecast with 2 policies as reference output
output_vars = [:forecastobs, :forecastpseudo]
modal_params = map(x -> x.value, m.parameters)
# sys_2pol = compute_system(m; tvis = true)
out1 = DSGE.forecast_one_draw(m, :mode, :full, output_vars, modal_params, df,
                              regime_switching = true, n_regimes = get_setting(m, :n_regimes))

# Test multiple alternative policies with multi-period altpolicies and temporary policies

## First start with fake temporary policy
temp_taylor_regime_eqcond_info = deepcopy(get_setting(m, :regime_eqcond_info))
temp_default_regime_eqcond_info = deepcopy(get_setting(m, :regime_eqcond_info))
for k in keys(temp_taylor_regime_eqcond_info)
    get_setting(m, :regime_eqcond_info)[k].weights =
        [get_setting(m, :regime_eqcond_info)[k].weights[1], get_setting(m, :regime_eqcond_info)[k].weights[2] / 2.,
         get_setting(m, :regime_eqcond_info)[k].weights[2] / 2.] # update the weights vector length
    temp_taylor_regime_eqcond_info[k] = DSGE.EqcondEntry(taylor_rule())
    temp_default_regime_eqcond_info[k] = DSGE.EqcondEntry(default_policy())
end
temp_taylor = MultiPeriodAltPolicy(:temporary_taylor, get_setting(m, :n_regimes), temp_taylor_regime_eqcond_info, gensys2 = true,
                                   temporary_altpolicy_names = [:taylor_rule],
                                   temporary_altpolicy_length = get_setting(m, :temporary_altpol_length),
                                   infoset = copy(get_setting(m, :tvis_information_set))) # also test w/ tvis and w/out
temp_default = MultiPeriodAltPolicy(:temporary_default, get_setting(m, :n_regimes), temp_default_regime_eqcond_info,
                                    gensys2 = false, temporary_altpolicy_names = [:default]) # also check without gensys2 on

delete!(DSGE.get_settings(m), :alternative_policies) # delete to update typing
m <= Setting(:alternative_policies, [DSGE.taylor_rule(), temp_taylor])
# sys_taylor = compute_system(m; tvis = true)
out_taylor_temp_taylor = DSGE.forecast_one_draw(m, :mode, :full, output_vars, modal_params, df,
                                                regime_switching = true, n_regimes = get_setting(m, :n_regimes))
m <= Setting(:alternative_policies, [DSGE.default_policy(), temp_default])
# sys_default = compute_system(m; tvis = true)
out_default_temp_default = DSGE.forecast_one_draw(m, :mode, :full, output_vars, modal_params, df,
                                                  regime_switching = true, n_regimes = get_setting(m, :n_regimes))

## Use implemented policy as nontrivial third policy
temp_flexait_zlb_regime_eqcond_info = deepcopy(get_setting(m, :regime_eqcond_info))
for k in keys(temp_taylor_regime_eqcond_info)
    get_setting(m, :regime_eqcond_info)[k].weights =
        [get_setting(m, :regime_eqcond_info)[k].weights[1] / 2., 2 * get_setting(m, :regime_eqcond_info)[k].weights[2],
         get_setting(m, :regime_eqcond_info)[k].weights[1] / 2.] # need to multiply by 2 b/c of how weights are set previously
    if get_setting(m, :regime_eqcond_info)[k].alternative_policy.key == :zero_rate
        temp_flexait_zlb_regime_eqcond_info[k] = EqcondEntry(zero_rate())
    else
        temp_flexait_zlb_regime_eqcond_info[k] = EqcondEntry(flexible_ait())
    end
end
temp_flexait_zlb = MultiPeriodAltPolicy(:temporary_flexait_zlb, get_setting(m, :n_regimes), temp_flexait_zlb_regime_eqcond_info, gensys2 = true,
                                        temporary_altpolicy_names = [:zero_rate],
                                        temporary_altpolicy_length = get_setting(m, :temporary_altpol_length),
                                        infoset = copy(get_setting(m, :tvis_information_set)))
m <= Setting(:alternative_policies, [DSGE.default_policy(), temp_flexait_zlb])
# sys_flexait = compute_system(m; tvis = true)
out_flexait_zlb_temp_flexait_zlb = DSGE.forecast_one_draw(m, :mode, :full, output_vars, modal_params, df,
                                                          regime_switching = true, n_regimes = get_setting(m, :n_regimes))

## Now add nontrivial temporary policy that is different from implemented policy
temp_flexible_ait_regime_eqcond_info = deepcopy(get_setting(m, :regime_eqcond_info))
for k in keys(temp_taylor_regime_eqcond_info)
    if get_setting(m, :regime_eqcond_info)[k].alternative_policy.key == :zero_rate
        temp_flexible_ait_regime_eqcond_info[k] = DSGE.EqcondEntry(flexible_ait()) # temporary flexible ait during ZLB periods
    else
        temp_flexible_ait_regime_eqcond_info[k] = DSGE.EqcondEntry(taylor_rule()) # followed by Taylor rule in the end
    end
end
temp_flexible_ait = MultiPeriodAltPolicy(:temporary_flexible_ait, get_setting(m, :n_regimes), temp_flexible_ait_regime_eqcond_info, gensys2 = true,
                                         temporary_altpolicy_names = [:flexible_ait],
                                         temporary_altpolicy_length = get_setting(m, :temporary_altpol_length),
                                         infoset = copy(get_setting(m, :tvis_information_set))) # also test w/ tvis and w/out

m <= Setting(:alternative_policies, [DSGE.default_policy(), temp_flexible_ait])
out_temp_flexible_ait = DSGE.forecast_one_draw(m, :mode, :full, output_vars, modal_params, df,
                                               regime_switching = true, n_regimes = get_setting(m, :n_regimes))

if regenerate_reference_forecasts
    h5open(joinpath(dirname(@__FILE__), "../reference/multiple_altpol_imperfect_awareness_output.h5"), "w") do file
        write(file, "forecastobs", out_temp_flexible_ait[:forecastobs])
        write(file, "forecastpseudo", out_temp_flexible_ait[:forecastpseudo])
    end
end


@testset "Multiple alternative policies, incl. MultiPeriodAltPolicy, with imperfect awareness" begin

    for k in keys(out_taylor_temp_taylor) # check implementations that should yield same results as a forecast
        @test out_taylor_temp_taylor[k] ≈ out1[k] # with 2 alternative policies and imperfect awareness
        @test out_default_temp_default[k] ≈ out1[k]
        @test out_flexait_zlb_temp_flexait_zlb[k] ≈ out1[k]
        @test out_temp_flexible_ait[k] ≈
            h5read(joinpath(dirname(@__FILE__), "../reference/multiple_altpol_imperfect_awareness_output.h5"), string(k))
    end

end
nothing
