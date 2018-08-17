using DSGE, Base.Test, HDF5, JLD

# Initialize model object
m = AnSchorfheide(testing = true)
m <= Setting(:forecast_horizons, 12)

estroot = normpath(joinpath(dirname(@__FILE__), "..", "reference"))
overrides = forecast_input_file_overrides(m)
overrides[:mode] = joinpath(estroot, "optimize.h5")

# Zero out measurement error (or else filtering shocks won't work exactly)
for para in [:e_y, :e_π, :e_R]
    m[para].value = 0
end

# Simulate dummy scenarios
def = Scenario(:defscen, "Test Default Scenario", [:obs_gdp, :obs_cpi], Symbol[], "REF")
alt = Scenario(:altscen, "Test Alternative Scenario", [:obs_gdp, :obs_cpi], [:g_sh, :rm_sh], "REF")

forecast_scenario(m, def, verbose = :none)
forecast_scenario(m, alt, verbose = :none)

expect = jldopen(get_scenario_input_file(m, alt), "r") do file
    read(file, "arr")
end
actual = jldopen(get_scenario_output_files(m, alt, [:forecastobs])[:forecastobs], "r") do file
    dict = read(file, "observable_indices")
    inds = map(var -> dict[var], def.target_names)
    read(file, "arr")[:, inds, :]
end
@testset "Test dummy scenarios" begin
    @test expect ≈ actual
end

# Simulate always switching
probs_enter = vcat([1], zeros(11))
probs_exit  = zeros(12)
allalt = SwitchingScenario(:allalt, alt, def, probs_enter, probs_exit)
simulate_switching(m, allalt, verbose = :none)

original_draws = h5read(get_scenario_output_files(m, alt, [:forecastobs])[:forecastobs], "arr")
actual = h5read(get_scenario_output_files(m, allalt, [:forecastobs])[:forecastobs], "arr")
@testset "Test switching scenarios" begin
    @test original_draws ≈ actual
end

# Simulate never switching
probs_enter = zeros(12)
probs_exit  = zeros(12)
alldef = SwitchingScenario(:alldef, alt, def, probs_enter, probs_exit)
simulate_switching(m, alldef, verbose = :none)

default_draws = h5read(get_scenario_output_files(m, def, [:forecastobs])[:forecastobs], "arr")
actual = h5read(get_scenario_output_files(m, alldef, [:forecastobs])[:forecastobs], "arr")
@testset "Test never switching scenarios" begin
    for i = 1:10 # Check that the ith draw of actual matches some draw j from default_draws
        @test for j = 1:10
            if actual[i, :, :] ≈ default_draws[j, :, :]
                return true
            end
        end
    end
end

# Simulate switching in at t = 5 and out at t = 9
probs_enter = vcat(zeros(4), [1], zeros(7))
probs_exit  = vcat(zeros(8), [1], zeros(3))
somealt = SwitchingScenario(:somealt, alt, def, probs_enter, probs_exit)
simulate_switching(m, somealt, verbose = :none)

@testset "Test switching scenarios, where the switch happens at particular points in time" begin
    actual = h5read(get_scenario_output_files(m, somealt, [:forecastobs])[:forecastobs], "arr")
    for i = 1:10
        @test actual[i, :, 5:8] ≈ original_draws[i, :, 1:4]
        @test for j = 1:10
            if actual[i, :, [1:4;9:12]] ≈ default_draws[j, :, [1:4;9:12]]
                return true
            end
        end
    end
end

# Transform single scenarios
scenario_means_bands(m, alt, verbose = :none)
scenario_means_bands(m, allalt, verbose = :none)
actual = squeeze(mean(original_draws, 1), 1)

@testset "Test single scenarios" begin
    mb1 = read_scenario_mb(m, alt, :forecastutobs)
    mb2 = read_scenario_mb(m, allalt, :forecastutobs)
    for mb in [mb1, mb2]
        for (var, ind) in mb.metadata[:indices]
            @test mb.means[var] ≈ actual[ind, :]
        end
    end

    mb1 = read_scenario_mb(m, alt, :forecastobs)
    mb2 = read_scenario_mb(m, allalt, :forecastobs)
    for mb in [mb1, mb2]
        @test mb.means[:obs_gdp] ≈ quartertoannual(actual[1, :])
        @test mb.means[:obs_cpi] ≈ quartertoannual(actual[2, :])
        @test mb.means[:obs_nominalrate] ≈ actual[3, :]
    end

    mb1 = read_scenario_mb(m, alt, :forecast4qobs)
    mb2 = read_scenario_mb(m, allalt, :forecast4qobs)
    for mb in [mb1, mb2]
        @test mb.means[:obs_gdp] ≈ loggrowthtopct_4q_approx(actual[1, :], zeros(3))
        @test mb.means[:obs_cpi] ≈ loggrowthtopct_4q_approx(actual[2, :], zeros(3))
        @test mb.means[:obs_nominalrate] ≈ actual[3, :]
    end
end

# Aggregate scenarios equally
aggall = ScenarioAggregate(:aggall, "Test Scenario Aggregate", AbstractScenario[def, alt],
                           [0.5, 0.5], 20, false, "REF")
scenario_means_bands(m, aggall, verbose = :none)
actual = squeeze(mean(cat(1, default_draws, original_draws), 1), 1)

@testset "Test scenario aggregation" begin
    mb = read_scenario_mb(m, aggall, :forecastutobs)
    for (var, ind) in mb.metadata[:indices]
        @test mb.means[var] ≈ actual[ind, :]
    end

    # Aggregate scenarios to only have alternative scenaro
    aggalt = ScenarioAggregate(:aggalt, "Test Scenario Aggregate", AbstractScenario[def, alt],
                               [0.0, 1.0], 20, false, "REF")
    scenario_means_bands(m, aggalt, verbose = :none)
    actual = squeeze(mean(repeat(original_draws, outer = [2, 1, 1]), 1), 1)

    mb = read_scenario_mb(m, aggalt, :forecastutobs)
    for (var, ind) in mb.metadata[:indices]
        @test mb.means[var] ≈ actual[ind, :]
    end
end

# Aggregate scenarios, drawing with replacement
aggrep1 = ScenarioAggregate(:aggrep1, "Test Scenario Aggregate", AbstractScenario[def, alt],
                            [0.5, 0.5], 20, true, "REF")
aggrep2 = ScenarioAggregate(:aggrep2, "Test Scenario Aggregate", AbstractScenario[def, alt],
                            [0.0, 1.0], 20, true, "REF")

# Simulate dummy scenario with shock scaling
scale = Scenario(:altscen, "Test Shock Scaling Scenario", [:obs_gdp, :obs_cpi], [:g_sh, :rm_sh], "REF",
                 shock_scaling = 2.0)
forecast_scenario(m, scale, verbose = :none)

expect = jldopen(get_scenario_input_file(m, scale), "r") do file
    2.0 * read(file, "arr")
end
actual = jldopen(get_scenario_output_files(m, alt, [:forecastobs])[:forecastobs], "r") do file
    dict = read(file, "observable_indices")
    inds = map(var -> dict[var], def.target_names)
    read(file, "arr")[:, inds, :]
end
@testset "Test scenario aggregation, drawing with replacement" begin
    @test expect ≈ actual
end

nothing
