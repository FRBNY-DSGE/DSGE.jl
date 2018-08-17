using DSGE, Base.Test

# Initialize model object
m = Model990(testing = true)
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m <= Setting(:use_population_forecast, true)

# Forecast under historical rule
params = [α.value for α in m.parameters]
df = load_data(m, verbose = :none)

output_vars = [:histpseudo, :forecastpseudo, :irfpseudo]
out_hist = DSGE.forecast_one_draw(m, :mode, :none, output_vars, params, df, verbose = :none)

# Assign alternative policy
m <= Setting(:alternative_policy, DSGE.taylor93())

# Compare equilibrium conditions under historical and alternative rules
Γ0_hist, Γ1_hist, _ = eqcond(m)
Γ0_alt,  Γ1_alt,  _ = alternative_policy(m).eqcond(m)

eq_mp    = m.equilibrium_conditions[:eq_mp]
eq_other = setdiff(1:n_equilibrium_conditions(m), eq_mp)

@testset "Compare AltPolicy eqconds under hist and alt rules" begin
    @test Γ0_hist[eq_other, :] ≈ Γ0_alt[eq_other, :]
    @test Γ1_hist[eq_other, :] ≈ Γ1_alt[eq_other, :]
    @test !(Γ0_hist[eq_mp, :]  ≈ Γ0_alt[eq_mp, :])
    @test !(Γ1_hist[eq_mp, :]  ≈ Γ1_alt[eq_mp, :])
end

# Test error thrown if trying to run shockdec products
@testset "Check AltPolicy error thrown with shockdec products" begin
    @test_throws ErrorException forecast_one(m, :mode, :none, [:shockdecobs])
    @test_throws ErrorException forecast_one(m, :mode, :none, [:dettrendobs])
    @test_throws ErrorException forecast_one(m, :mode, :none, [:trendobs])
end

# Forecast under Taylor 93
out_alt = DSGE.forecast_one_draw(m, :mode, :none, output_vars, params, df, verbose = :none)

@testset "Check forecast under Taylor 93" begin
    @test out_hist[:histpseudo] ≈ out_alt[:histpseudo]
    @test !(out_hist[:forecastpseudo] ≈ out_alt[:forecastpseudo])
end
