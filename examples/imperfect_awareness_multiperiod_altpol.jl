using DSGE, ModelConstructors, CSV, DataFrames, Dates, Test

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

# Initialize model objects
m = Model1002("ss62")
m <= Setting(:date_forecast_start,  fcast_date)
m <= Setting(:date_conditional_end, fcast_date)
m <= Setting(:forecast_horizons, subtract_quarters(date_fcast_end, fcast_date))

# Load data
df = CSV.read(joinpath(dirname(@__FILE__), "..", "test", "reference", "ss62_modal_forecast_data.csv"), DataFrame)

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

# Define a function to help update the mappings
# from model regimes to parameter regimes
function set_regime_vals_fnct(m, n)
    m2p = get_setting(m, :model2para_regime)

    for (k, v) in m2p
        for i in 7:n
            v[i] = 1 # For regimes 7 onward, use parameter regime 1
        end
    end

    m
end

# Set up credibility for ZLB and alternative policy
m <= Setting(:gensys2, true)
m <= Setting(:replace_eqcond, true)
m <= Setting(:uncertain_temporary_altpolicy, true)
m <= Setting(:uncertain_altpolicy, true)

# Find the first regime in which gensys2 should apply
gensys2_first_regime = findfirst([get_setting(m, :regime_dates)[i] .== start_zlb_date for i in 1:get_setting(m, :n_regimes)])
n_zlb_reg = DSGE.subtract_quarters(end_zlb_date, start_zlb_date) + 1
if end_zlb_date < get_setting(m, :regime_dates)[get_setting(m, :n_regimes)]
    set_regime_vals_fnct(m, get_setting(m, :n_regimes))
else
    # technically the second input is `gensys2_first_regime + (n_zlb_reg - 1) + 1` b/c
    # gensys2_first_regime is the first regime for the temporary ZLB, so we need to subtract 1
    # from n_zlb_reg, but we need to also add 1 for the lift-off regime
    set_regime_vals_fnct(m, gensys2_first_regime + n_zlb_reg)
end

# Construct regime_eqcond_info
reg_dates = deepcopy(get_setting(m, :regime_dates))
regime_eqcond_info = Dict{Int, EqcondEntry}()
weights = [end_tvcred_level, 1. - end_tvcred_level] # dummy values for now, time-varying credibility added later
for (regind, date) in zip(gensys2_first_regime:(n_zlb_reg - 1 + gensys2_first_regime), # See comments starting at line 75
                          DSGE.quarter_range(reg_dates[gensys2_first_regime], # looping over the date range of the ZLB
                                                iterate_quarters(reg_dates[gensys2_first_regime], n_zlb_reg - 1)))
    reg_dates[regind] = date # updating the regime date dictionary
    regime_eqcond_info[regind] = EqcondEntry(zlb_rule(), weights) # specify which regimes have a ZLB rule applying
end
reg_dates[n_zlb_reg + gensys2_first_regime] = iterate_quarters(reg_dates[gensys2_first_regime], n_zlb_reg) # add liftoff regime date
regime_eqcond_info[n_zlb_reg + gensys2_first_regime] = EqcondEntry(flexible_ait(), weights) # liftoff policy is flexible AIT
m <= Setting(:regime_dates,               reg_dates)
m <= Setting(:regime_eqcond_info,         regime_eqcond_info)
m <= Setting(:temporary_altpolicy_length, n_zlb_reg) # length of the temporary ZLB, otherwise a guess is made. Since we have time-varying credibility, this guess would be wrong
m <= Setting(:zlb_rule_value,             .1) # (optional) set the interest rate during the ZLB, defaults to forecast_zlb_value(m) if not available
m <= Setting(:temporary_altpolicy_names,  [:zlb_rule]) # Ensures compute_system knows the regimes w/ a zlb_rule in EqcondEntry are gensys2 regimes

# Set up regime indices
setup_regime_switching_inds!(m; cond_type = cond_type)

# Set up time-varying credibility (which extends beyond the ZLB, so we need to add more regimes)
tvcred_dates = (start_zlb_date, DSGE.iterate_quarters(start_zlb_date, tvcred_n_qtrs))
n_tvcred_reg = DSGE.subtract_quarters(tvcred_dates[2], tvcred_dates[1]) + 1 # number of regimes for time-varying credibility
reg_start = DSGE.subtract_quarters(tvcred_dates[1], start_zlb_date) + gensys2_first_regime # first regime for time-varying credibility

# Add extra entries for regime_eqcond_info with dummy credibility weights
for (regind, date) in zip(reg_start:(reg_start + n_tvcred_reg),
                          DSGE.quarter_range(reg_dates[reg_start], DSGE.iterate_quarters(reg_dates[reg_start], n_tvcred_reg - 1)))
    reg_dates[regind] = date
    if date > end_zlb_date # Add extra regimes with flexible AIT as the policy after the ZLB's end b/c credibility is still changing
        regime_eqcond_info[regind] = DSGE.EqcondEntry(flexible_ait(), [end_tvcred_level, 1. - end_tvcred_level])
    end
end
m <= Setting(:regime_dates,             reg_dates)
m <= Setting(:regime_eqcond_info, regime_eqcond_info)
setup_regime_switching_inds!(m; cond_type = cond_type)

# Update mapping of model regimes to parameter regimes
set_regime_vals_fnct(m, get_setting(m, :n_regimes))

# Set up TVIS information set
m <= Setting(:tvis_information_set, vcat([i:i for i in 1:(gensys2_first_regime - 1)],
                                         [i:17 for i in gensys2_first_regime:17],
                                         [i:i for i in 18:get_setting(m, :n_regimes)]))

# Now actually populate the weights in regime_eqcond_info so that
# credibility/awareness increases linearly over time
credvec = collect(range(start_tvcred_level, stop = end_tvcred_level, length = n_tvcred_reg))
for (i, k) in enumerate(sort!(collect(keys(regime_eqcond_info))))
    if tvcred_dates[1] <= reg_dates[k] <= tvcred_dates[2]
        get_setting(m, :regime_eqcond_info)[k].weights = [credvec[i], 1. - credvec[i]]
    end
end

# We implement a shorter ZLB with a return to the historical
# monetary policy reaction function as the other policy
# agents believe could occur. This is where the `MultiPeriodAltPolicy` comes into play.
altpol_eqcond_info = deepcopy(get_setting(m, :regime_eqcond_info))
for (reg, eq_entry) in altpol_eqcond_info
    if reg >= start_reg_baseline + 5
        altpol_eqcond_info[reg] = DSGE.EqcondEntry(DSGE.default_policy(), [1.])
    end
    altpol_eqcond_info[reg].weights = [1.]
end

# The key for the MultiPeriodAltPolicy can be anything; it is never used in the code
# to infer some property. Of course, it should still be a descriptive name.
altpol = MultiPeriodAltPolicy(:other_altpol, get_setting(m, :n_regimes), altpol_eqcond_info, gensys2 = true,
                              temporary_altpolicy_names = [:zlb_rule],
                              temporary_altpolicy_length = 6,
                              infoset = copy(get_setting(m, :tvis_information_set)))
delete!(m.settings, :alternative_policies) # if this setting already exists, the type may be wrong
m <= Setting(:alternative_policies, DSGE.AbstractAltPolicy[altpol])

# Run forecasts
output_vars = [:histobs, :histpseudo, :forecastpseudo, :forecastobs]
modal_params = map(x -> x.value, m.parameters)
fcast = DSGE.forecast_one_draw(m, :mode, cond_type, output_vars, modal_params, df,
                               regime_switching = true, n_regimes = get_setting(m, :n_regimes))

nothing
