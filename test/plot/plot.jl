using Plots
path = dirname(@__FILE__)

# Initialize the plotting backend
gr()

# Initialize model object
m = AnSchorfheide(testing = true)
m <= Setting(:saveroot, tempdir())
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m <= Setting(:date_conditional_end, quartertodate("2015-Q4"))
m <= Setting(:use_population_forecast, true)
m <= Setting(:forecast_horizons, 12)
m <= Setting(:impulse_response_horizons, 8)

# Run full-distribution forecast
estroot = normpath(joinpath(dirname(@__FILE__), "..", "reference"))
overrides = forecast_input_file_overrides(m)
overrides[:mode] = joinpath(estroot, "optimize.h5")
overrides[:full] = joinpath(estroot, "metropolis_hastings.h5")
output_vars = add_requisite_output_vars([:histobs, :forecastobs, :shockdecobs, :irfobs])

@everywhere using DSGE
m <= Setting(:forecast_block_size, 5)
forecast_one(m, :full, :none, output_vars, verbose = :none)
compute_meansbands(m, :full, :none, output_vars; verbose = :none)
println("The following warning is expected test behavior:")

# Plot history and forecast
plot_history_and_forecast(m, :obs_nominalrate, :obs, :full, :none,
                          bdd_and_unbdd = true,
                          start_date = DSGE.quartertodate("2007-Q1"),
                          verbose = :none)

# Plot forecast comparison
plot_forecast_comparison(m, m, :obs_nominalrate, :obs, :full, :none,
                         bdd_and_unbdd = true,
                         start_date = DSGE.quartertodate("2007-Q1"),
                         verbose = :none)

# Plot shock decomposition
plot_shock_decomposition(m, :obs_nominalrate, :obs, :full, :none, verbose = :none)

# Plot IRF
plot_impulse_response(m, :rm_sh, collect(keys(m.observables)), :obs, :full, :none,
                      verbose = :none)
plot_impulse_response(m, :rm_sh, collect(keys(m.observables)), :obs, :full, :none,
                      flip = true, verbose = :none)
plot_prior_posterior(m, :τ)
plot_prior_posterior(m, [:τ, :σ_g])
cloud = load("$path/../reference/output_data/an_schorfheide/ss0/estimate/raw/smc_cloud_vint=151001.jld2")["cloud"]
plot_posterior_intervals(m, cloud = cloud)
plot_posterior_intervals(m, [cloud, cloud])

# Plot scenario, zeroing out measurement error
for para in [:e_y, :e_π, :e_R]
    m[para].value = 0
end
alt = Scenario(:altscen, "Test Alternative Scenario", [:obs_gdp, :obs_cpi], [:g_sh, :rm_sh], "REF")
forecast_scenario(m, alt, verbose = :none)
scenario_means_bands(m, alt, verbose = :none)
plot_scenario(m, :obs_nominalrate, :obs, alt, untrans = true, verbose = :none)
plot_scenario(m, :obs_nominalrate, :obs, alt, verbose = :none)
plot_scenario(m, :obs_nominalrate, :obs, alt, fourquarter = true, verbose = :none)
@testset "Test scenario exception" begin
    @test_throws ErrorException plot_scenario(m, :obs_nominalrate, :obs, alt,
                                              untrans = true, fourquarter = true)
end

# Hair plot
realized = load_data(m, verbose = :none)
hist_mb = read_mb(m, :full, :none, :histobs)
fcast_mb = read_mb(m, :full, :none, :bddforecastobs)
hair_plot(:obs_nominalrate, realized, [hist_mb], [fcast_mb];
          plotroot = saveroot(m), verbose = :none)

# Things requiring estimated files
tmp_saveroot = saveroot(m)
m <= Setting(:saveroot, "$path/../reference")
plot_history_and_forecast(m, :y_t, :pseudo, :full, :none,
                          bdd_and_unbdd = true,
                          start_date = DSGE.quartertodate("2007-Q1"),
                          verbose = :none, plotroot = tmp_saveroot)
plot_forecast_comparison(m, m, :y_t, :pseudo, :full, :none,
                         bdd_and_unbdd = true,
                         start_date = DSGE.quartertodate("2007-Q1"),
                         verbose = :none, plotroot = tmp_saveroot)
plot_shock_decomposition(m, :y_t, :pseudo, :full, :none, verbose = :none, plotroot = tmp_saveroot)
plot_impulse_response(m, :rm_sh, collect(keys(m.pseudo_observables)), :pseudo, :full, :none,
                      verbose = :none, plotroot = tmp_saveroot)
plot_impulse_response(m, :rm_sh, collect(keys(m.pseudo_observables)), :pseudo, :full, :none,
                      flip = true, verbose = :none, plotroot = tmp_saveroot)
plot_altpolicies([m], :obs_gdp, :obs, :none, plotroot = tmp_saveroot)
plot_altpolicies([m], :y_t, :pseudo, :none, plotroot = tmp_saveroot)
plot_forecast_comparison(m, m, :obs_gdp, :obs, :full, :none, plotroot = tmp_saveroot)
plot_forecast_comparison(m, m, :y_t, :pseudo, :full, :none, plotroot = tmp_saveroot)
plot_forecast_decomposition(m, m, :obs_gdp, :obs, :full, :none, :none, plotroot = tmp_saveroot)
plot_forecast_decomposition(m, m, :y_t, :pseudo, :full, :none, :none, plotroot = tmp_saveroot)
plot_forecast_sequence([m,m], [:mode, :mode], [:none, :none], m, :mode, :none, :obs, :obs_gdp,
                       plotroot = tmp_saveroot, start_date = DSGE.quartertodate("1960-Q1"),
                       end_date = DSGE.quartertodate("1961-Q1"))
plot_forecast_sequence([m,m], [:full, :full], [:none, :none], m, :full, :none, :obs, :obs_gdp,
                       plotroot = tmp_saveroot, start_date = DSGE.quartertodate("1960-Q1"),
                       end_date = DSGE.quartertodate("1961-Q1"))
plot_forecast_sequence([m,m], [:mode, :mode], [:none, :none], m, :mode, :none, :pseudo, :y_t,
                       plotroot = tmp_saveroot, start_date = DSGE.quartertodate("1960-Q1"),
                       end_date = DSGE.quartertodate("1961-Q1"))
plot_forecast_sequence([m,m], [:full, :full], [:none, :none], m, :full, :none, :pseudo, :y_t,
                       plotroot = tmp_saveroot, start_date = DSGE.quartertodate("1960-Q1"),
                       end_date = DSGE.quartertodate("1961-Q1"))
