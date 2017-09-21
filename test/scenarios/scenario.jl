using DSGE, Base.Test

# Initialize model object
m = AnSchorfheide(testing = true)
m <= Setting(:forecast_horizons, 2)

# Scenario
scen = Scenario(:testscen, "Test Scenario", [:obs_gdp, :obs_cpi], [:g_sh, :rm_sh], "REF")
@test scen.shock_scaling == 1.0
@test n_targets(scen) == 2
@test n_instruments(scen) == 2

for var in scen.target_names
    scen.targets[var] = rand(2)
end
@test n_target_horizons(scen) == 2

act = targets_to_data(m, scen)
exp = copy(scen.targets)
exp[:date] = DSGE.quarter_range(date_forecast_start(m), date_forecast_end(m))
exp[:obs_nominalrate] = fill(NaN, 2)
for col in names(exp)
    @test act[col] == exp[col] || (all(isnan(act[col])) && all(isnan(exp[col])))
end

# Scenario with shock scaling
scale = Scenario(:scaledscen, "Test Shock Scaling Scenario", [:obs_gdp, :obs_cpi], [:g_sh, :rm_sh], "REF",
                 shock_scaling = 2.0)
@test scale.shock_scaling == 2.0

# SwitchingScenario
@test_throws ErrorException SwitchingScenario(:switch, scen, scen, [-1.0, 1.0], [1.0, 1.0])
@test_throws ErrorException SwitchingScenario(:switch, scen, scen, [1.0, 1.0], [-1.0, 1.0])
@test_throws ErrorException SwitchingScenario(:switch, scen, scen, [1.0], [1.0, 1.0])
@test_throws ErrorException SwitchingScenario(:switch, scen, scen, ones(3), zeros(3))

# ScenarioAggregate
@test_throws ErrorException ScenarioAggregate(:agg, "Test Scenario Aggregate", [SingleScenario[scen]],
                                              [0.5, 0.5], 10, true, "REF")
@test_throws ErrorException ScenarioAggregate(:agg, "Test Scenario Aggregate", [SingleScenario[scen]],
                                              [-1.0], 10, true, "REF")
@test_throws ErrorException ScenarioAggregate(:agg, "Test Scenario Aggregate", [SingleScenario[scen]],
                                              [0.5], 10, true, "REF")