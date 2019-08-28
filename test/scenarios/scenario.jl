using DSGE, Test

# Initialize model object
m = AnSchorfheide(testing = true)
m <= Setting(:forecast_horizons, 2)

# Scenario
scen = Scenario(:testscen, "Test Scenario", [:obs_gdp, :obs_cpi], [:g_sh, :rm_sh], "REF")
@testset "Test scenario basic constructor" begin
    @test scen.shock_scaling == 1.0
    @test n_targets(scen) == 2
    @test n_instruments(scen) == 2

    for var in scen.target_names
        scen.targets[!, var] = rand(2)
    end
    @test n_target_horizons(scen) == 2

    global act = targets_to_data(m, scen)
    global expt = copy(scen.targets)
    expt[!, :date] = DSGE.quarter_range(date_forecast_start(m), date_forecast_end(m))
    expt[!, :obs_nominalrate] = fill(missing, 2)
    for col in names(expt)
        global check = false
        if ismissing(act[!, col] == expt[!, col])
            check = (all(ismissing.(act[!, col])) && all(ismissing.(expt[!, col])))
        else
            check = (act[!, col] == expt[!, col])
        end
        @test check
    end
end

# Scenario with shock scaling
@testset "Test scenario constructor with shock scaling" begin
    scale = Scenario(:scaledscen, "Test Shock Scaling Scenario", [:obs_gdp, :obs_cpi], [:g_sh, :rm_sh], "REF",
                     shock_scaling = 2.0)
    @test scale.shock_scaling == 2.0

    # SwitchingScenario
    @test_throws ErrorException SwitchingScenario(:switch, scen, scen, [1.0], [1.0, 1.0])
    @test_throws ErrorException SwitchingScenario(:switch, scen, scen, [-1.0, 1.0], [1.0, 1.0])
    @test_throws ErrorException SwitchingScenario(:switch, scen, scen, [1.0, 1.0], [-1.0, 1.0])

    # ScenarioAggregate
    @test_throws ErrorException ScenarioAggregate(:agg, "Test Scenario Aggregate", AbstractScenario[scen],
                                                  [0.5, 0.5], 10, true, "REF")
    @test_throws ErrorException ScenarioAggregate(:agg, "Test Scenario Aggregate", AbstractScenario[scen],
                                                  [-1.0], 10, true, "REF")
    @test_throws ErrorException ScenarioAggregate(:agg, "Test Scenario Aggregate", AbstractScenario[scen],
                                                  [0.5], 10, true, "REF")
end
