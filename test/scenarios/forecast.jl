using DSGE, Test

writing_output = false

# Initialize model object
m = AnSchorfheide(testing = true)
m <= Setting(:forecast_horizons, 12)

# Define dummy scenario
scen = Scenario(:testscen, "Test Scenario", [:obs_gdp, :obs_cpi], [:g_sh, :rm_sh], "REF")

# Check compute_scenario_system throws error before zeroing out measurement error
@testset "Check compute_scenario_system error throwing" begin
    @test_throws ErrorException compute_scenario_system(m, scen)
end

# Zero out measurement error (or else filtering shocks won't work exactly)
for para in [:e_y, :e_π, :e_R]
    m[para].value = 0
end

# Test compute_scenario_system
@testset "Test compute_scenario_system, filter_shocks!, and shock_scaling" begin
    global sys = compute_scenario_system(m, scen)
    @test all(x -> x == 0, sys[:CCC])
    @test all(x -> x == 0, sys[:DD])
    @test all(x -> x == 0, sys[:DD_pseudo])
    for shock in keys(m.exogenous_shocks)
        ind = m.exogenous_shocks[shock]
        if shock in scen.instrument_names
            @test sys[:QQ][ind, ind] != 0
        else
            @test sys[:QQ][ind, ind] == 0
        end
    end

    # Test filter_shocks!
    for var in scen.target_names
        scen.targets[!, var] = rand(forecast_horizons(m))
    end
    global forecastshocks = filter_shocks!(m, scen, sys)
    for var in scen.instrument_names
        ind = m.exogenous_shocks[var]
        @test scen.instruments[!, var] == forecastshocks[ind, :]
    end

    s_T = zeros(n_states_augmented(m))
    global _, forecastobs,_ ,_  = forecast(m, sys, s_T, shocks = forecastshocks)
    for var in scen.target_names
        ind = m.observables[var]
        @test scen.targets[!, var] ≈ forecastobs[ind, :]
    end

    # Test scenario with shock scaling
    global scale_scenario = Scenario(:scaledscen, "Test Shock Scaling Scenario",
                            [:obs_gdp, :obs_cpi], [:g_sh, :rm_sh], "REF",
                            shock_scaling = 2.0)
    scale_scenario.targets = scen.targets

    global forecastshocks = filter_shocks!(m, scale_scenario, sys)
    for var in scale_scenario.target_names
        ind = m.observables[var]
        @test scale_scenario.targets[!, var] == scen.targets[!, var]
    end
    for var in scale_scenario.instrument_names
        ind = m.exogenous_shocks[var]
        @test scale_scenario.instruments[!, var] == forecastshocks[ind, :]
    end
end

###########################
# forecast_scenario
###########################
#=
file = "$path/../reference/forecast_scenario_draw.jld2"
sys = compute_scenario_system(m, scen)

forecast_output = DSGE.forecast_scenario_draw(m, scen, sys, 1)

if writing_output
    JLD2.jldopen(file, true, true, true, IOStream) do file
        file["forecast_output"] = forecast_output
    end
end
forecast_output_test = load(file, "forecast_output")

@testset "Test forecast_scenario_draw, forecast_scenario" begin
    @test forecast_output == forecast_output_test
end
=#
