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

# Set up credibility for ZLB and alternative policy, if applicable
m <= Setting(:uncertain_temporary_altpolicy, true)
# m <= Setting(:temporary_altpolicy, :zero_rate)
m <= Setting(:uncertain_altpolicy, true)
gensys2_first_regime = get_setting(m, :n_regimes) + 1
get_setting(m, :regime_dates)[gensys2_first_regime] = start_zlb_date
n_zlb_reg = DSGE.subtract_quarters(end_zlb_date, start_zlb_date) + 1

# Set up regime_eqcond_info
m <= Setting(:replace_eqcond, true)
reg_dates = deepcopy(get_setting(m, :regime_dates))
regime_eqcond_info = Dict{Int, DSGE.EqcondEntry}()
for (regind, date) in zip(gensys2_first_regime:(n_zlb_reg - 1 + gensys2_first_regime), # See comments starting at line 57
                          DSGE.quarter_range(reg_dates[gensys2_first_regime],
                                             DSGE.iterate_quarters(reg_dates[gensys2_first_regime], n_zlb_reg - 1)))
    reg_dates[regind] = date
    regime_eqcond_info[regind] = DSGE.EqcondEntry(DSGE.zero_rate(), [0., 1.])
end
reg_dates[n_zlb_reg + gensys2_first_regime] = DSGE.iterate_quarters(reg_dates[gensys2_first_regime], n_zlb_reg)
regime_eqcond_info[n_zlb_reg + gensys2_first_regime] = DSGE.EqcondEntry(DSGE.flexible_ait(), [0., 1.])
nreg0 = length(reg_dates)
m <= Setting(:regime_dates,             reg_dates)
m <= Setting(:regime_eqcond_info, regime_eqcond_info)

# Set up regime indices and extra regimes for parameters
setup_regime_switching_inds!(m; cond_type = :full)
set_regime_vals_fnct(m, get_setting(m, :n_regimes))

# Set up TVIS information set
m <= Setting(:tvis_information_set, vcat([i:i for i in 1:(gensys2_first_regime - 1)],
                                         [i:get_setting(m, :n_regimes) for i in
                                          gensys2_first_regime:get_setting(m, :n_regimes)]))

# output_vars = [:histobs, :histpseudo, :forecastobs, :forecastpseudo, :forecaststates]
output_vars = [:forecastobs, :forecastpseudo]
modal_params = map(x -> x.value, m.parameters)
@assert !haskey(m.settings, :temporary_altpolicy_length)
outp0 = DSGE.forecast_one_draw(m, :mode, :full, output_vars, modal_params, df,
                                     regime_switching = true, n_regimes = get_setting(m, :n_regimes))
for i in keys(get_setting(m, :regime_eqcond_info))
    get_setting(m, :regime_eqcond_info)[i].weights = [θ[:cred], 1. - θ[:cred]]
end

# Compare fixed vs. time-varying credibility
outp33 = DSGE.forecast_one_draw(m, :mode, :full, output_vars, modal_params, df,
                                regime_switching = true, n_regimes = get_setting(m, :n_regimes))
sysp33 = compute_system(m; tvis = true)
for i in nreg0:(nreg0 + 15)
    reg_dates[i] = DSGE.iterate_quarters(reg_dates[nreg0], i - nreg0)
    regime_eqcond_info[i] = DSGE.EqcondEntry(DSGE.flexible_ait(), [.33, 1-.33])
end
m <= Setting(:regime_dates,             reg_dates)
m <= Setting(:regime_eqcond_info, regime_eqcond_info)
setup_regime_switching_inds!(m; cond_type = :full)
m <= Setting(:tvis_information_set, vcat([i:i for i in 1:(gensys2_first_regime - 1)],
                                         [i:get_setting(m, :n_regimes) for i in
                                          gensys2_first_regime:get_setting(m, :n_regimes)]))
set_regime_vals_fnct(m, get_setting(m, :n_regimes))
m <= Setting(:temporary_altpolicy_length, n_zlb_reg)
for i in keys(get_setting(m, :regime_eqcond_info))
    get_setting(m, :regime_eqcond_info)[i].weights = [.33, 1-.33]
end
sysp33_tv = compute_system(m; tvis = true)
outp33_tv = DSGE.forecast_one_draw(m, :mode, :full, output_vars, modal_params, df,
                                   regime_switching = true, n_regimes = get_setting(m, :n_regimes))

@testset "Compare Fixed to Time-Varying Credibility" begin
    for reg in 1:n_regimes(sysp33)
        @test sysp33[reg, :TTT] ≈ sysp33_tv[reg, :TTT]
        @test sysp33[reg, :RRR] ≈ sysp33_tv[reg, :RRR]
        @test sysp33[reg, :CCC] ≈ sysp33_tv[reg, :CCC]
        @test sysp33[reg, :QQ] ≈ sysp33_tv[reg, :QQ]
        @test sysp33[reg, :ZZ] ≈ sysp33_tv[reg, :ZZ]
        @test sysp33[reg, :DD] ≈ sysp33_tv[reg, :DD]
    end
    for k in keys(outp33)
        @test outp33[k] ≈ outp33_tv[k]
    end

    if !regenerate_reference_forecasts
        testfcast = h5read(joinpath(dirname(@__FILE__), "..", "reference", "tvcred_reference_meas_fix_forecast.h5"), "forecastobs")
        inds = vcat(1:9, 12:13, 20:21)  # ignore k-periods ahead observables, reference data generated when a bug existed
        @test maximum(abs.(testfcast[inds, :] -
                           outp33[:forecastobs][inds, :])) < 1.1e-3#5e-4 # some numerical differences b/c fixes to calculations
    end

end

# Compare perfect credibility to permanent case
for i in keys(get_setting(m, :regime_eqcond_info))
    get_setting(m, :regime_eqcond_info)[i].weights = [1., 0.]
end
credvec = ones(17)
for (i, k) in enumerate(sort!(collect(keys(regime_eqcond_info))))
    if (get_setting(m, :temporary_altpolicy_length) - 1) < i
        get_setting(m, :regime_eqcond_info)[k].weights = [credvec[i - (get_setting(m, :temporary_altpolicy_length) - 1)],
                                                          1. - credvec[i - (get_setting(m, :temporary_altpolicy_length) - 1)]]
    end
end
out_temp = DSGE.forecast_one_draw(m, :mode, :full, output_vars, modal_params, df,
                                  regime_switching = true, n_regimes = get_setting(m, :n_regimes))

m <= Setting(:uncertain_temporary_altpolicy, false)
m <= Setting(:uncertain_altpolicy, false)

credvec = collect(range(0., stop = 1., length = 17))
for (i, k) in enumerate(sort!(collect(keys(regime_eqcond_info))))
    if i <= length(credvec)
        get_setting(m, :regime_eqcond_info)[k].weights = [credvec[i], 1.0 - credvec[i]]
    end
end

out_perm = DSGE.forecast_one_draw(m, :mode, :full, output_vars, modal_params, df,
                              regime_switching = true, n_regimes = get_setting(m, :n_regimes))

# Check perfectly credible ZLB match permanent equivalent
m <= Setting(:uncertain_temporary_altpolicy, false)
m <= Setting(:uncertain_altpolicy, true)
credvec = collect(range(0., stop = 1., length = 17))
for (i, k) in enumerate(sort!(collect(keys(regime_eqcond_info))))
    if (get_setting(m, :temporary_altpolicy_length) ) < i
        get_setting(m, :regime_eqcond_info)[k].weights = [credvec[i - (get_setting(m, :temporary_altpolicy_length))],
                                                          1. - credvec[i - (get_setting(m, :temporary_altpolicy_length))]]
    end
end

out_credzlb = DSGE.forecast_one_draw(m, :mode, :full, output_vars, modal_params, df,
                                     regime_switching = true, n_regimes = get_setting(m, :n_regimes))

m <= Setting(:uncertain_altpolicy, true)
m <= Setting(:uncertain_temporary_altpolicy, true)
for i in keys(get_setting(m, :regime_eqcond_info))
    get_setting(m, :regime_eqcond_info)[i].weights = [1., 0.]
end
credvec = collect(range(0., stop = 1., length = 17))

for (i, k) in enumerate(sort!(collect(keys(regime_eqcond_info))))
    if (get_setting(m, :temporary_altpolicy_length)) < i
        get_setting(m, :regime_eqcond_info)[k].weights = [credvec[i - (get_setting(m, :temporary_altpolicy_length))],
                                                          1. - credvec[i - (get_setting(m, :temporary_altpolicy_length))]]
    else
        firsted = 1.0
        get_setting(m, :regime_eqcond_info)[k].weights = [firsted, 1.0 - firsted]
    end
end

out_temp_credzlb = DSGE.forecast_one_draw(m, :mode, :full, output_vars, modal_params, df,
                                          regime_switching = true, n_regimes = get_setting(m, :n_regimes))


if !regenerate_reference_forecasts
    @testset "Compare Perfect Credibility to Permanent Alternative Policy" begin
        @test out_temp[:forecastobs] ≈ out_perm[:forecastobs]
        @test out_temp_credzlb[:forecastobs] ≈ out_credzlb[:forecastobs]
    end
end

for i in keys(get_setting(m, :regime_eqcond_info))
    get_setting(m, :regime_eqcond_info)[i].weights = [0., 1.]
end
credvec = collect(range(0., stop = 1., length = 17))
for (i, k) in enumerate(sort!(collect(keys(regime_eqcond_info))))
    if (get_setting(m, :temporary_altpolicy_length) - 1) < i
        get_setting(m, :regime_eqcond_info)[k].weights = [credvec[i - (get_setting(m, :temporary_altpolicy_length) - 1)],
                                                          1. - credvec[i - (get_setting(m, :temporary_altpolicy_length) - 1)]]
    end
end
out1 = DSGE.forecast_one_draw(m, :mode, :full, output_vars, modal_params, df,
                              regime_switching = true, n_regimes = get_setting(m, :n_regimes))

if !regenerate_reference_forecasts
    inds = vcat(1:9, 12:13)
    @testset "Compare TV Credibility to Reference Forecast" begin
        tvtestfcast = h5read(joinpath(dirname(@__FILE__), "..", "reference", "tvcred_reference_forecast.h5"),
                             VERSION >= v"1.5" ? "tvforecastobs_1p5" : "tvforecastobs")
        @test maximum(abs.(out1[:forecastobs][inds, :] - tvtestfcast[inds, :])) < 1e-3
    end
end



if regenerate_reference_forecasts
    h5open(joinpath(dirname(@__FILE__), "..", "reference", "tvcred_reference_meas_fix_forecast.h5"), "w") do file
        write(file, "para", θ10)
        write(file, "forecastobs", outp33[:forecastobs])
        if VERSION >= v"1.5"
            write(file, "tvforecastobs", tvtestfcast_othervers)
            write(file, "tvforecastobs_1p5", out1[:forecastobs])
        else
            write(file, "forecastobs_tv0to1", out1[:forecastobs])
            write(file, "forecastobs_fixed33", outp33[:forecastobs])
            write(file, "forecastobs_fixed0", outp0[:forecastobs])
            write(file, "forecastobs_fixed100", out_temp[:forecastobs])
            write(file, "forecastobs_tv0to1_ZLBcred_1", out_credzlb[:forecastobs])
        end
    end
end

nothing
