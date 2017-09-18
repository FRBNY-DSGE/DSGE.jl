using DSGE, Base.Test

# Initialize model object
m = AnSchorfheide(testing = true)
m <= Setting(:forecast_horizons, 12)

# Define dummy scenario
scen = Scenario(:testscen, "Test Scenario", [:obs_gdp, :obs_cpi], [:g_sh, :rm_sh], "REF")

# Check compute_scenario_system throws error before zeroing out measurement error
@test_throws ErrorException compute_scenario_system(m, scen)

# Zero out measurement error (or else filtering shocks won't work exactly)
for para in [:e_y, :e_π, :e_R]
    m[para].value = 0
end

# Test compute_scenario_system
sys = compute_scenario_system(m, scen)
@test all(x -> x == 0, sys[:CCC])
@test all(x -> x == 0, sys[:DD])
@test all(x -> x == 0, sys[:DD_pseudo])
for shock in keys(m.exogenous_shocks)
    i = m.exogenous_shocks[shock]
    if shock in scen.instrument_names
        @test sys[:QQ][i, i] != 0
    else
        @test sys[:QQ][i, i] == 0
    end
end

# Test filter_shocks!
for var in scen.target_names
    scen.targets[var] = rand(forecast_horizons(m))
end
forecastshocks = filter_shocks!(m, scen, sys)
for var in scen.instrument_names
    i = m.exogenous_shocks[var]
    @test scen.instruments[var] == forecastshocks[i, :]
end

s_T = zeros(n_states_augmented(m))
_, forecastobs, _, _ = forecast(m, sys, s_T, shocks = forecastshocks)
for var in scen.target_names
    i = m.observables[var]
    @test scen.targets[var] ≈ forecastobs[i, :]
end