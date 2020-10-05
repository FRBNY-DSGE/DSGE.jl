# Test that regime-switch in 2020-Q3 to flexible AIT rule runs properly

using Test, ModelConstructors, DSGE, Dates, FileIO, Random

generate_output_data = false

m30 = Model1002("ss30")
m30 <= Setting(:date_forecast_start, Date(2020, 9, 30))
m30 <= Setting(:date_conditional_end, Date(2020, 9, 30))
df = load_data(m30; check_empty_columns = false)
sys30 = compute_system(m30)

# Compute as is


# Now check m10 generates same results when using flexible AIT policy as an alternative policy
m10 = Model1002("ss10"; custom_settings = Dict{Symbol, Setting}(:add_pgap => Setting(:add_pgap, true),
                                                                :add_ygap => Setting(:add_ygap, true)))
m10 <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m10), 2 => Date(2020, 9, 30)))
m10 <= Setting(:regime_switching, true)
DSGE.setup_regime_switching_inds!(m10)
m10 <= Setting(:alternative_policy, DSGE.flexible_ait())
m10 <= Setting(:pgap_type, :flexible_ait)
m10 <= Setting(:pgap_value, 0.)
m10 <= Setting(:ygap_type, :flexible_ait)
m10 <= Setting(:ygap_value, 12.)
m10 <= Setting(:ait_Thalf, 10.)
m10 <= Setting(:gdp_Thalf, 10.)
m10 <= Setting(:flexible_ait_ρ_smooth, 0.)
m10 <= Setting(:flexible_ait_φ_π, 6.)
m10 <= Setting(:flexible_ait_φ_y, 6.)
sys10_noalt = compute_system(m10; apply_altpolicy = false)
sys10_alt = compute_system(m10; apply_altpolicy = true)

sys30[1, :TTT] ≈ sys10_noalt[1, :TTT]
!(sys30[2, :TTT] ≈ sys10_noalt[2, :TTT])
sys30[1, :TTT] ≈ sys10_alt[1, :TTT]
sys30[2, :TTT] ≈ sys10_alt[2, :TTT]
