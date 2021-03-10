using DSGE, ModelConstructors, JLD2, CSV, DataFrames, Dates, Test

regenerate_reference_output = false

# Data vintages and settings
fcast_date      = DSGE.quartertodate("2021-Q1")
date_fcast_end  = iterate_quarters(fcast_date, 40)
cond_type       = :full

# Flexible AIT Rule
tvcred_n_qtrs          = 24 # length of tv cred period
start_tvcred_level     = 0.
end_tvcred_level       = 1.
φ_π                    = 15.
φ_y                    = 3.
Thalf                  = 10.
ρ_smooth               = .75
pgap_init              = 0.125 # NOTE: positive pgap and ygap values imply
ygap_init              = 12.   # NEGATIVE initial pgap and ygaps (i.e. ygap_init = 12 implies
pgap_ygap_init_date    = Date(2020, 6, 30)  # an initial ygap of -12 in 2020:Q2)
start_zlb_date         = Date(2020, 12, 31) # start and end dates should be inclusive, i.e. first period w/zero rate
end_zlb_date           = Date(2023, 9, 30)  # is start_zlb_date, and the last period is end_zlb_date

# additional settings to implement Flexible AIT rule
flexait_custom = Dict{Symbol, Setting}(:forecast_horizons => Setting(:forecast_horizons,
                                                                     subtract_quarters(date_fcast_end, fcast_date)),
                                       :add_iid_cond_obs_gdp_meas_err =>
                                       Setting(:add_iid_cond_obs_gdp_meas_err, true),
                                       :add_iid_cond_obs_corepce_meas_err =>
                                       Setting(:add_iid_cond_obs_corepce_meas_err, true),
                                       :add_pseudo_corepce => Setting(:add_pseudo_corepce, true))

# Initialize model objects
m = Model1002("ss62", custom_settings = flexait_custom)
m <= Setting(:date_forecast_start,  fcast_date)
m <= Setting(:date_conditional_end, fcast_date)

m2p = get_setting(m, :model2para_regime)
m2p[:σ_condgdp] = Dict(1 => 1, 2 => 1, 3 => 1, 4 => 1, 5 => 1, 6 => 2)
m2p[:σ_condcorepce] = Dict(1 => 1, 2 => 1, 3 => 1, 4 => 1, 5 => 1, 6 => 2)
set_regime_val!(m[:σ_condgdp], 1, 0.)
set_regime_val!(m[:σ_condgdp], 2, 2. * m[:σ_gdp].value)
set_regime_val!(m[:σ_condcorepce], 1, 0.)
set_regime_val!(m[:σ_condcorepce], 2, 1.25 * m[:σ_corepce].value)

# Load data
df = CSV.read(joinpath(dirname(@__FILE__), "..", "reference", "ss62_modal_forecast_data.csv"), DataFrame)

# Set up Flexible AIT with a temporary ZLB. Note that df will be altered
start_reg_baseline = 6

θ = (φ_π = φ_π, φ_y = φ_y, Thalf = Thalf, pgap = pgap_init,
     ygap = ygap_init, ρ_smooth = ρ_smooth,
     historical_policy = default_policy()) # NamedTuple with parameters of the Flexible AIT policy
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

function set_regime_vals_fnct(m, n)
    m2p = get_setting(m, :model2para_regime)

    for (k, v) in m2p
        for i in 7:n
            v[i] = 1 # For regimes 7 onward, use parameter regime 1
        end
    end

    m
end

m <= Setting(:gensys2, true)

# Set up credibility for ZLB and alternative policy, if applicable
m <= Setting(:uncertain_temporary_altpolicy, true)
m <= Setting(:uncertain_altpolicy, true)

# Find the first regime in which gensys2 should apply
gensys2_first_regime = findfirst([get_setting(m, :regime_dates)[i] .== start_zlb_date for i in 1:get_setting(m, :n_regimes)])
n_zlb_reg = DSGE.subtract_quarters(end_zlb_date, start_zlb_date) + 1
if end_zlb_date < get_setting(m, :regime_dates)[get_setting(m, :n_regimes)]
    set_regime_vals_fnct(m, get_setting(m, :n_regimes))
else
    # technically gensys2_first_regime + (n_zlb_reg - 1) + 1 b/c
    # gensys2_first_regime is the first regime for the temporary ZLB, so need to subtract 1
    # from n_zlb_reg, but need to also add 1 for the lift-off regime
    set_regime_vals_fnct(m, gensys2_first_regime + n_zlb_reg)
end

m <= Setting(:replace_eqcond, true)
reg_dates = deepcopy(get_setting(m, :regime_dates))
regime_eqcond_info = Dict{Int, DSGE.EqcondEntry}()
weights = [end_tvcred_level, 1. - end_tvcred_level]
for (regind, date) in zip(gensys2_first_regime:(n_zlb_reg - 1 + gensys2_first_regime), # See comments starting at line 57
                          DSGE.quarter_range(reg_dates[gensys2_first_regime],
                                                iterate_quarters(reg_dates[gensys2_first_regime], n_zlb_reg - 1)))
    reg_dates[regind] = date
    regime_eqcond_info[regind] = DSGE.EqcondEntry(DSGE.zero_rate(), weights)
end
reg_dates[n_zlb_reg + gensys2_first_regime] = iterate_quarters(reg_dates[gensys2_first_regime], n_zlb_reg)
regime_eqcond_info[n_zlb_reg + gensys2_first_regime] = DSGE.EqcondEntry(flexible_ait(), weights)
m <= Setting(:regime_dates,               reg_dates)
m <= Setting(:regime_eqcond_info,         regime_eqcond_info)
m <= Setting(:temporary_altpolicy_length, n_zlb_reg)

# Remove anticipated MP shocks upon entering ZLB
m <= Setting(:remove_rm_shocks, gensys2_first_regime)

# Set up regime indices
setup_regime_switching_inds!(m; cond_type = cond_type)

# Set up TVIS information set
m <= Setting(:tvis_information_set, vcat([i:i for i in 1:(gensys2_first_regime - 1)],
                                         [i:get_setting(m, :n_regimes) for i in
                                          gensys2_first_regime:get_setting(m, :n_regimes)]))

tvcred_dates = (start_zlb_date, DSGE.iterate_quarters(start_zlb_date, tvcred_n_qtrs))
n_tvcred_reg = DSGE.subtract_quarters(tvcred_dates[2], tvcred_dates[1]) + 1 # number of regimes for time-varying credibility
reg_start = DSGE.subtract_quarters(tvcred_dates[1], start_zlb_date) + gensys2_first_regime # first regime for time-varying credibility
for (regind, date) in zip(reg_start:(reg_start + n_tvcred_reg),
                          DSGE.quarter_range(reg_dates[reg_start], DSGE.iterate_quarters(reg_dates[reg_start], n_tvcred_reg - 1)))
    reg_dates[regind] = date
    if date > end_zlb_date
        regime_eqcond_info[regind] = DSGE.EqcondEntry(flexible_ait(), [end_tvcred_level, 1. - end_tvcred_level])
    end
end
m <= Setting(:regime_dates,             reg_dates)
m <= Setting(:regime_eqcond_info, regime_eqcond_info)
setup_regime_switching_inds!(m; cond_type = cond_type)
set_regime_vals_fnct(m, get_setting(m, :n_regimes))
m <= Setting(:tvis_information_set, vcat([i:i for i in 1:(gensys2_first_regime - 1)],
                                         [i:get_setting(m, :n_regimes) for i in
                                          gensys2_first_regime:get_setting(m, :n_regimes)]))
credvec = collect(range(start_tvcred_level, stop = end_tvcred_level, length = n_tvcred_reg))
for (i, k) in enumerate(sort!(collect(keys(regime_eqcond_info))))
    if tvcred_dates[1] <= reg_dates[k] <= tvcred_dates[2]
        get_setting(m, :regime_eqcond_info)[k].weights = [credvec[i], 1. - credvec[i]]
    end
end

m <= Setting(:zero_rate_zlb_value, 0.1)

# Run forecasts
output_vars = [:histobs, :histpseudo, :forecastpseudo, :forecastobs]
modal_params = map(x -> x.value, m.parameters)

# Create two-rule as alternative policy
tworule_eqcond_info = deepcopy(get_setting(m, :regime_eqcond_info))
for (reg, eq_entry) in tworule_eqcond_info
    if reg >= start_reg_baseline + 5
        tworule_eqcond_info[reg] = DSGE.EqcondEntry(DSGE.default_policy(), [1.])
    end
    tworule_eqcond_info[reg].weights = [1.]
end
tworule = MultiPeriodAltPolicy(:two_rule, get_setting(m, :n_regimes), tworule_eqcond_info, gensys2 = true,
                               temporary_altpolicy_names = [:zero_rate],
                               temporary_altpolicy_length = 6,
                               infoset = copy(get_setting(m, :tvis_information_set)))
delete!(DSGE.get_settings(m), :alternative_policies)
m <= Setting(:alternative_policies, DSGE.AbstractAltPolicy[tworule])

sys1 = compute_system(m; tvis = true)
using BenchmarkTools
@btime begin
    sys1 = compute_system(m; tvis = true)
    nothing
end
m <= Setting(:perfect_cred_regime_mapping, Dict(i => 18 for i in 18:29))
tworule2 = MultiPeriodAltPolicy(:two_rule, get_setting(m, :n_regimes), tworule_eqcond_info, gensys2 = true,
                                temporary_altpolicy_names = [:zero_rate],
                                temporary_altpolicy_length = 6,
                                perfect_cred_regime_mapping = Dict(i => 11 for i in 11:29),
                                infoset = copy(get_setting(m, :tvis_information_set)))
m <= Setting(:alternative_policies, DSGE.AbstractAltPolicy[tworule2])
sys2 = compute_system(m; tvis = true)
@btime begin
    sys2 = compute_system(m; tvis = true)
    nothing
end
@assert false
fcast = DSGE.forecast_one_draw(m, :mode, cond_type, output_vars, modal_params, df,
                               regime_switching = true, n_regimes = get_setting(m, :n_regimes))

if regenerate_reference_output
    JLD2.jldopen(joinpath(dirname(@__FILE__), "..", "reference", "ss62_modal_forecast_output.jld2"), true, true, true, IOStream) do file
        write(file, "forecastobs", fcast[:forecastobs])
        write(file, "forecastpseudo", fcast[:forecastpseudo])
        write(file, "histpseudo", fcast[:histpseudo])
    end
end

@testset "Compare Modal Forecast for Model 1002 ss62 against reference" begin
    ref_out = JLD2.jldopen(joinpath(dirname(@__FILE__), "..", "reference", "ss62_modal_forecast_output.jld2"), "r")
    @test fcast[:forecastobs] ≈ ref_out["forecastobs"]
    @test fcast[:forecastpseudo] ≈ ref_out["forecastpseudo"]
    @test fcast[:histpseudo] ≈ ref_out["histpseudo"]
end

nothing
