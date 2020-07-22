using DSGE, Test, ModelConstructors

# Initialize model object
m = Model1002(testing = true)
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m <= Setting(:use_population_forecast, true)
m <= Setting(:data_vintage, "160812")

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

# Repeat for alt_inflation
m <= Setting(:alternative_policy, DSGE.alt_inflation())
Γ0_alt_inf,  Γ1_alt_inf,  ~ = alternative_policy(m).eqcond(m)

# Repeat for rw
m <= Setting(:alternative_policy, DSGE.rw())
Γ0_alt_rw,  Γ1_alt_rw,  ~ = alternative_policy(m).eqcond(m)

# Repeat for rw_zero_rate
m <= Setting(:alternative_policy, DSGE.rw_zero_rate())
Γ0_alt_rw_zero,  Γ1_alt_rw_zero,  ~ = alternative_policy(m).eqcond(m)

# Repeat for zero_rate
m <= Setting(:alternative_policy, DSGE.zero_rate())
Γ0_alt_zero,  Γ1_alt_zero,  ~ = alternative_policy(m).eqcond(m)

# Repeat for smooth_ait_gdp
m <= Setting(:alternative_policy, DSGE.smooth_ait_gdp())
Γ0_alt_smooth,  Γ1_alt_smooth,  ~ = alternative_policy(m).eqcond(m)

# Repeat for smooth_ait_gdp_alt
m <= Setting(:alternative_policy, DSGE.smooth_ait_gdp_alt())
Γ0_alt_smooth_alt,  Γ1_alt_smooth_alt,  ~ = alternative_policy(m).eqcond(m)

# Repeat for ait
m <= Setting(:alternative_policy, DSGE.ait())
m <= Setting(:pgap_type, :ait)
m <= Setting(:pgap_value, 0.)
Γ0_alt_ait,  Γ1_alt_ait,  ~ = alternative_policy(m).eqcond(m)

# Repeat for ngdp
m <= Setting(:alternative_policy, DSGE.ngdp())
m <= Setting(:pgap_type, :ngdp)
m <= Setting(:pgap_value, 12.)
Γ0_alt_ngdp,  Γ1_alt_ngdp,  ~ = alternative_policy(m).eqcond(m)



# Getting indices
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

    @test Γ0_hist[eq_other, :] ≈ Γ0_alt_inf[eq_other, :]
    @test Γ1_hist[eq_other, :] ≈ Γ1_alt_inf[eq_other, :]
    @test !(Γ0_hist[eq_mp, :]  ≈ Γ0_alt_inf[eq_mp, :])
    @test !(Γ1_hist[eq_mp, :]  ≈ Γ1_alt_inf[eq_mp, :])

    @test Γ0_hist[eq_other, :] ≈ Γ0_alt_rw[eq_other, :]
    @test Γ1_hist[eq_other, :] ≈ Γ1_alt_rw[eq_other, :]
    @test !(Γ0_hist[eq_mp, :]  ≈ Γ0_alt_rw[eq_mp, :])
    @test !(Γ1_hist[eq_mp, :]  ≈ Γ1_alt_rw[eq_mp, :])

    @test Γ0_hist[eq_other, :] ≈ Γ0_alt_rw_zero[eq_other, :]
    @test Γ1_hist[eq_other, :] ≈ Γ1_alt_rw_zero[eq_other, :]
    @test !(Γ0_hist[eq_mp, :]  ≈ Γ0_alt_rw_zero[eq_mp, :])
    @test !(Γ1_hist[eq_mp, :]  ≈ Γ1_alt_rw_zero[eq_mp, :])

    @test Γ0_hist[eq_other, :] ≈ Γ0_alt_zero[eq_other, :]
    @test Γ1_hist[eq_other, :] ≈ Γ1_alt_zero[eq_other, :]
    @test !(Γ0_hist[eq_mp, :]  ≈ Γ0_alt_zero[eq_mp, :])
    @test !(Γ1_hist[eq_mp, :]  ≈ Γ1_alt_zero[eq_mp, :])

    @test Γ0_hist[eq_other, :] ≈ Γ0_alt_smooth[eq_other, :]
    @test Γ1_hist[eq_other, :] ≈ Γ1_alt_smooth[eq_other, :]
    @test !(Γ0_hist[eq_mp, :]  ≈ Γ0_alt_smooth[eq_mp, :])
    @test !(Γ1_hist[eq_mp, :]  ≈ Γ1_alt_smooth[eq_mp, :])

    @test Γ0_hist[eq_other, :] ≈ Γ0_alt_smooth_alt[eq_other, :]
    @test Γ1_hist[eq_other, :] ≈ Γ1_alt_smooth_alt[eq_other, :]
    @test !(Γ0_hist[eq_mp, :]  ≈ Γ0_alt_smooth_alt[eq_mp, :])
    @test !(Γ1_hist[eq_mp, :]  ≈ Γ1_alt_smooth_alt[eq_mp, :])

    @test Γ0_hist[eq_other, :] ≈ Γ0_alt_ait[eq_other, :]
    @test Γ1_hist[eq_other, :] ≈ Γ1_alt_ait[eq_other, :]
    @test !(Γ0_hist[eq_mp, :]  ≈ Γ0_alt_ait[eq_mp, :])
    @test !(Γ1_hist[eq_mp, :]  ≈ Γ1_alt_ait[eq_mp, :])

    @test Γ0_hist[eq_other, :] ≈ Γ0_alt_ngdp[eq_other, :]
    @test Γ1_hist[eq_other, :] ≈ Γ1_alt_ngdp[eq_other, :]
    @test !(Γ0_hist[eq_mp, :]  ≈ Γ0_alt_ngdp[eq_mp, :])
    @test !(Γ1_hist[eq_mp, :]  ≈ Γ1_alt_ngdp[eq_mp, :])
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

# Forecast under alt_inflation and rw
m <= Setting(:alternative_policy, DSGE.alt_inflation())
out_alt_inf = DSGE.forecast_one_draw(m, :mode, :none, output_vars, paras, df, verbose = :none)
m <= Setting(:alternative_policy, DSGE.rw())
out_alt_rw = DSGE.forecast_one_draw(m, :mode, :none, output_vars, paras, df, verbose = :none)

m <= Setting(:alternative_policy, DSGE.rw_zero_rate())
out_alt_rw_zero = DSGE.forecast_one_draw(m, :mode, :none, output_vars, paras, df, verbose = :none)
m <= Setting(:alternative_policy, DSGE.zero_rate())
out_alt_zero = DSGE.forecast_one_draw(m, :mode, :none, output_vars, paras, df, verbose = :none)

m <= Setting(:alternative_policy, DSGE.smooth_ait_gdp())
out_alt_smooth = DSGE.forecast_one_draw(m, :mode, :none, output_vars, paras, df, verbose = :none)
m <= Setting(:alternative_policy, DSGE.smooth_ait_gdp_alt())
out_alt_smooth_alt = DSGE.forecast_one_draw(m, :mode, :none, output_vars, paras, df, verbose = :none)

m <= Setting(:alternative_policy, DSGE.ait())
out_alt_ait = DSGE.forecast_one_draw(m, :mode, :none, output_vars, paras, df, verbose = :none)
m <= Setting(:alternative_policy, DSGE.ngdp())
out_alt_ngdp = DSGE.forecast_one_draw(m, :mode, :none, output_vars, paras, df, verbose = :none)

@testset "Check forecast under Taylor 93" begin
    @test out_hist[:histpseudo] ≈ out_alt93[:histpseudo]
    @test !(out_hist[:forecastpseudo] ≈ out_alt93[:forecastpseudo])
end
@testset "Check forecast under Taylor 99" begin
    @test out_hist[:histpseudo] ≈ out_alt99[:histpseudo]
    @test !(out_hist[:forecastpseudo] ≈ out_alt99[:forecastpseudo])
end

@testset "Check forecast under alt_inflation" begin
    @test out_hist[:histpseudo] ≈ out_alt_inf[:histpseudo]
    @test !(out_hist[:forecastpseudo] ≈ out_alt_inf[:forecastpseudo])
end
@testset "Check forecast under rw" begin
    @test out_hist[:histpseudo] ≈ out_alt_rw[:histpseudo]
    @test !(out_hist[:forecastpseudo] ≈ out_alt_rw[:forecastpseudo])
end

@testset "Check forecast under rw_zero_rate" begin
    @test out_hist[:histpseudo] ≈ out_alt_rw_zero[:histpseudo]
    @test !(out_hist[:forecastpseudo] ≈ out_alt_rw_zero[:forecastpseudo])
end
@testset "Check forecast under zero_rate" begin
    @test out_hist[:histpseudo] ≈ out_alt_zero[:histpseudo]
    @test !(out_hist[:forecastpseudo] ≈ out_alt_zero[:forecastpseudo])
end

@testset "Check forecast under smooth_ait_gdp" begin
    @test out_hist[:histpseudo] ≈ out_alt_smooth[:histpseudo]
    @test !(out_hist[:forecastpseudo] ≈ out_alt_smooth[:forecastpseudo])
end
@testset "Check forecast under smooth_ait_gdp_alt" begin
    @test out_hist[:histpseudo] ≈ out_alt_smooth_alt[:histpseudo]
    @test !(out_hist[:forecastpseudo] ≈ out_alt_smooth_alt[:forecastpseudo])
end

@testset "Check forecast under ait" begin
    @test out_hist[:histpseudo] ≈ out_alt_ait[:histpseudo]
    @test !(out_hist[:forecastpseudo] ≈ out_alt_ait[:forecastpseudo])
end
@testset "Check forecast under ngdp" begin
    @test out_hist[:histpseudo] ≈ out_alt_ngdp[:histpseudo]
    @test !(out_hist[:forecastpseudo] ≈ out_alt_ngdp[:forecastpseudo])
end


# Checking if setting variables works
## These values are random and not actually good choices
m <= Setting(:ait_Thalf, 14.2)
m <= Setting(:gdp_Thalf, 5.6)
m <= Setting(:smooth_ait_gdp_ρ_smooth, 0.23)
m <= Setting(:smooth_ait_gdp_φ_π, 0.36)
m <= Setting(:smooth_ait_gdp_φ_y, 0.56)
m <= Setting(:smooth_ait_gdp_alt_ρ_smooth, 0.76)
m <= Setting(:smooth_ait_gdp_alt_φ_π, 0.11)
m <= Setting(:smooth_ait_gdp_alt_φ_y, 0.89)
m <= Setting(:rw_ρ_smooth, 0.45)
m <= Setting(:rw_φ_π, 0.83)
m <= Setting(:rw_φ_y, 0.94)
m <= Setting(:ρ_rw, 0.63)

@testset "Checking Settings assignments" begin
    @test get_setting(m, :ait_Thalf) == 14.2
    @test get_setting(m, :gdp_Thalf) == 5.6
    @test get_setting(m, :smooth_ait_gdp_ρ_smooth) == 0.23
    @test get_setting(m, :smooth_ait_gdp_φ_π) == 0.36
    @test get_setting(m, :smooth_ait_gdp_φ_y) == 0.56
    @test get_setting(m, :smooth_ait_gdp_alt_ρ_smooth) == 0.76
    @test get_setting(m, :smooth_ait_gdp_alt_φ_π) == 0.11
    @test get_setting(m, :smooth_ait_gdp_alt_φ_y) == 0.89
    @test get_setting(m, :rw_ρ_smooth) == 0.45
    @test get_setting(m, :rw_φ_π) == 0.83
    @test get_setting(m, :rw_φ_y) == 0.94
    @test get_setting(m, :ρ_rw) == 0.63
end
