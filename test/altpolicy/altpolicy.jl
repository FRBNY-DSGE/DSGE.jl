using DSGE, Test

# Initialize model object
m = Model990(testing = true)
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m <= Setting(:use_population_forecast, true)

# Forecast under historical rule
paras = [α.value for α in m.parameters]
df = load_data(m, verbose = :none)

output_vars = [:histpseudo, :forecastpseudo, :irfpseudo]
out_hist = DSGE.forecast_one_draw(m, :mode, :none, output_vars, paras, df, verbose = :none)

# Assign alternative policy
m <= Setting(:alternative_policy, DSGE.taylor93())

# Compare equilibrium conditions under historical and alternative rules
Γ0_hist, Γ1_hist, ~ = eqcond(m)
Γ0_alt93,  Γ1_alt93,  ~ = alternative_policy(m).eqcond(m)

# Repeat for Taylor99
m <= Setting(:alternative_policy, DSGE.taylor99())
Γ0_alt99,  Γ1_alt99,  ~ = alternative_policy(m).eqcond(m)

eq_mp    = m.equilibrium_conditions[:eq_mp]
eq_other = setdiff(1:n_equilibrium_conditions(m), eq_mp)

@testset "Compare AltPolicy eqconds under hist and alt rules" begin
    @test Γ0_hist[eq_other, :] ≈ Γ0_alt93[eq_other, :]
    @test Γ1_hist[eq_other, :] ≈ Γ1_alt93[eq_other, :]
    @test !(Γ0_hist[eq_mp, :]  ≈ Γ0_alt93[eq_mp, :])
    @test !(Γ1_hist[eq_mp, :]  ≈ Γ1_alt93[eq_mp, :])
    @test Γ0_hist[eq_other, :] ≈ Γ0_alt99[eq_other, :]
    @test Γ1_hist[eq_other, :] ≈ Γ1_alt99[eq_other, :]
    @test !(Γ0_hist[eq_mp, :]  ≈ Γ0_alt99[eq_mp, :])
    @test !(Γ1_hist[eq_mp, :]  ≈ Γ1_alt99[eq_mp, :])
end

# Test error thrown if trying to run shockdec products
@testset "Check AltPolicy error thrown with shockdec products" begin
    @test_throws ErrorException forecast_one(m, :mode, :none, [:shockdecobs])
    @test_throws ErrorException forecast_one(m, :mode, :none, [:dettrendobs])
    @test_throws ErrorException forecast_one(m, :mode, :none, [:trendobs])
end

# Forecast under Taylor 93 and Taylor 99
m <= Setting(:alternative_policy, DSGE.taylor93())
out_alt93 = DSGE.forecast_one_draw(m, :mode, :none, output_vars, paras, df, verbose = :none)
m <= Setting(:alternative_policy, DSGE.taylor99())
out_alt99 = DSGE.forecast_one_draw(m, :mode, :none, output_vars, paras, df, verbose = :none)

@testset "Check forecast under Taylor 93" begin
    @test out_hist[:histpseudo] ≈ out_alt93[:histpseudo]
    @test !(out_hist[:forecastpseudo] ≈ out_alt93[:forecastpseudo])
end
@testset "Check forecast under Taylor 99" begin
    @test out_hist[:histpseudo] ≈ out_alt99[:histpseudo]
    @test !(out_hist[:forecastpseudo] ≈ out_alt99[:forecastpseudo])
end
