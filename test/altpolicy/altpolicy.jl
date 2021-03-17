using DSGE, Test, ModelConstructors

# Initialize model object
m = Model990()
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m <= Setting(:use_population_forecast, true)
m990 = Model990(testing = true)
m990 <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m990 <= Setting(:use_population_forecast, true)

# Forecast under historical rule
paras = [α.value for α in m.parameters]
df = load_data(m990, verbose = :none)

output_vars = [:histpseudo, :forecastpseudo, :irfpseudo]
sys_hist = compute_system(m)
out_hist = DSGE.forecast_one_draw(m, :mode, :none, output_vars, paras, df, verbose = :none)

# Getting indices
eq_mp    = m.equilibrium_conditions[:eq_mp]
eq_other = setdiff(1:n_equilibrium_conditions(m), eq_mp)

# Assign alternative policy
setup_permanent_altpol!(m, DSGE.taylor93())

# Compare equilibrium conditions under historical and alternative rules
Γ0_hist, Γ1_hist, ~ = eqcond(m)
Γ0_alt93,  Γ1_alt93,  ~ = get_setting(m, :regime_eqcond_info)[2].alternative_policy.eqcond(m)

# Repeat for Taylor99
setup_permanent_altpol!(m, DSGE.taylor99())
Γ0_alt99,  Γ1_alt99,  ~ = get_setting(m, :regime_eqcond_info)[2].alternative_policy.eqcond(m)

## First, check if coefficients in matrices correct
### Taylor 93
endo = m.endogenous_states
@testset "Check Taylor93 rule has the right coefficients." begin
    @test Γ0_alt93[eq_mp, endo[:R_t]] == 1
    @test Γ0_alt93[eq_mp, endo[:π_a_t]] == -1.5/4
    @test Γ0_alt93[eq_mp, endo[:y_t]]  == -0.5/4
    @test Γ0_alt93[eq_mp, endo[:y_f_t]] == 0.5/4
end

### Taylor 99
@testset "Check Taylor99 rule has the right coefficients." begin
    @test Γ0_alt99[eq_mp, endo[:R_t]] == 1
    @test Γ0_alt99[eq_mp, endo[:π_a_t]] == -1.5/4
    @test Γ0_alt99[eq_mp, endo[:y_t]]  == -1/4
    @test Γ0_alt99[eq_mp, endo[:y_f_t]] == 1/4
end

# Testing
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

# Forecast under Taylor 93 and Taylor 99
setup_permanent_altpol!(m, DSGE.taylor93())
out_alt93 = DSGE.forecast_one_draw(m, :mode, :none, output_vars, paras, df, verbose = :none,
                                   regime_switching = true, n_regimes = get_setting(m, :n_regimes))
setup_permanent_altpol!(m, DSGE.taylor99())
out_alt99 = DSGE.forecast_one_draw(m, :mode, :none, output_vars, paras, df, verbose = :none,
                                   regime_switching = true, n_regimes = get_setting(m, :n_regimes))
m <= Setting(:date_forecast_start, DSGE.iterate_quarters(date_forecast_start(m), -1))
m <= Setting(:date_conditional_end, date_forecast_start(m))
setup_permanent_altpol!(m, DSGE.taylor99(); cond_type = :full)
out_alt99_cond = DSGE.forecast_one_draw(m, :mode, :full, output_vars, paras, df, verbose = :none,
                                        regime_switching = true, n_regimes = get_setting(m, :n_regimes))
m <= Setting(:date_forecast_start, DSGE.iterate_quarters(date_forecast_start(m), 1))
m <= Setting(:date_conditional_end, date_forecast_start(m))
@testset "Check forecast under Taylor 93" begin
    @test out_hist[:histpseudo] ≈ out_alt93[:histpseudo]
    @test !(out_hist[:forecastpseudo] ≈ out_alt93[:forecastpseudo])
end
@testset "Check forecast under Taylor 99" begin
    @test out_hist[:histpseudo] ≈ out_alt99[:histpseudo]
    @test !(out_hist[:forecastpseudo] ≈ out_alt99[:forecastpseudo])
end

####################
# New Policy Rules #
####################

m = Model1002("ss10"; custom_settings = Dict{Symbol, Setting}(:add_altpolicy_pgap => Setting(:add_altpolicy_pgap, true)))
my = Model1002("ss10"; custom_settings = Dict{Symbol, Setting}(:add_altpolicy_pgap => Setting(:add_altpolicy_pgap, true),
                                                               :add_altpolicy_ygap => Setting(:add_altpolicy_ygap, true)))
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m <= Setting(:use_population_forecast, true)
my <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
my <= Setting(:use_population_forecast, true)

# Getting indices
eq_mp    = m.equilibrium_conditions[:eq_mp]
eq_pgap  = m.equilibrium_conditions[:eq_pgap]
eq_other = setdiff(1:n_equilibrium_conditions(m), [eq_mp, eq_pgap])
eq_mpy    = my.equilibrium_conditions[:eq_mp]
eq_pgapy  = my.equilibrium_conditions[:eq_pgap]
eq_ygapy  = my.equilibrium_conditions[:eq_ygap]
eq_othery = setdiff(1:n_equilibrium_conditions(my), [eq_mpy, eq_pgapy, eq_ygapy])
eq_otherymp = 1:n_equilibrium_conditions(my)

# Compare equilibrium conditions under historical and alternative rules
Γ0_hist, Γ1_hist, ~ = eqcond(m)
Γ0_histy, Γ1_histy, ~ = eqcond(my)

# Matrices for alt_inflation
setup_permanent_altpol!(m, DSGE.alt_inflation())
Γ0_alt_inf,  Γ1_alt_inf,  ~ = get_setting(m, :regime_eqcond_info)[2].alternative_policy.eqcond(m)

# Repeat for zero_rate
setup_permanent_altpol!(m, DSGE.zero_rate())
Γ0_alt_zero,  Γ1_alt_zero,  ~ = get_setting(m, :regime_eqcond_info)[2].alternative_policy.eqcond(m)

# Repeat for smooth_ait_gdp
setup_permanent_altpol!(my, DSGE.smooth_ait_gdp())
Γ0_alt_smooth,  Γ1_alt_smooth,  ~ = get_setting(my, :regime_eqcond_info)[2].alternative_policy.eqcond(my)

# Repeat for smooth_ait_gdp_alt
setup_permanent_altpol!(my, DSGE.smooth_ait_gdp_alt())
Γ0_alt_smooth_alt,  Γ1_alt_smooth_alt,  ~ = get_setting(my, :regime_eqcond_info)[2].alternative_policy.eqcond(my)

# Repeat for flexible_ait
setup_permanent_altpol!(my, DSGE.flexible_ait())
Γ0_alt_flex_ait,  Γ1_alt_flex_ait,  ~ = get_setting(my, :regime_eqcond_info)[2].alternative_policy.eqcond(my)

# Repeat for ait
setup_permanent_altpol!(my, DSGE.ait())
my <= Setting(:pgap_type, :ait)
my <= Setting(:pgap_value, 0.)
Γ0_alt_ait,  Γ1_alt_ait,  ~ = get_setting(my, :regime_eqcond_info)[2].alternative_policy.eqcond(my)

# Repeat for ngdp
setup_permanent_altpol!(my, DSGE.ait())
my <= Setting(:pgap_type, :ngdp)
my <= Setting(:pgap_value, 12.)
Γ0_alt_ngdp,  Γ1_alt_ngdp,  ~ = get_setting(my, :regime_eqcond_info)[2].alternative_policy.eqcond(my)

# Repeat for rw
setup_permanent_altpol!(my, DSGE.rw())
Γ0_alt_rw,  Γ1_alt_rw,  ~ = get_setting(my, :regime_eqcond_info)[2].alternative_policy.eqcond(my)

# Repeat for rw_zero_rate
setup_permanent_altpol!(my, DSGE.rw_zero_rate())
Γ0_alt_rw_zero,  Γ1_alt_rw_zero,  ~ = get_setting(my, :regime_eqcond_info)[2].alternative_policy.eqcond(my)

# Testing
@testset "Compare AltPolicy eqconds under hist and alt rules" begin
    @test Γ0_hist[eq_other, :] ≈ Γ0_alt_inf[eq_other, :]
    @test Γ1_hist[eq_other, :] ≈ Γ1_alt_inf[eq_other, :]
    @test !(Γ0_hist[eq_mp, :]  ≈ Γ0_alt_inf[eq_mp, :])
    @test !(Γ1_hist[eq_mp, :]  ≈ Γ1_alt_inf[eq_mp, :])
    @test Γ0_hist[eq_other, eq_other] ≈ Γ0_alt_zero[eq_other, eq_other]
    @test Γ1_hist[eq_other, eq_other] ≈ Γ1_alt_zero[eq_other, eq_other]
    @test !(Γ0_hist[eq_mp, :]  ≈ Γ0_alt_zero[eq_mp, :])
    @test !(Γ1_hist[eq_mp, :]  ≈ Γ1_alt_zero[eq_mp, :])

    @test Γ0_histy[eq_othery, eq_otherymp] ≈ Γ0_alt_smooth[eq_othery, eq_otherymp]
    @test Γ1_histy[eq_othery, eq_otherymp] ≈ Γ1_alt_smooth[eq_othery, eq_otherymp]
    @test !(Γ0_histy[eq_mp, eq_other]  ≈ Γ0_alt_smooth[eq_mp, eq_othery])
    @test !(Γ1_histy[eq_mp, eq_other]  ≈ Γ1_alt_smooth[eq_mp, eq_othery])

    @test Γ0_histy[eq_othery, eq_otherymp] ≈ Γ0_alt_smooth_alt[eq_othery, eq_otherymp]
    @test Γ1_histy[eq_othery, eq_otherymp] ≈ Γ1_alt_smooth_alt[eq_othery, eq_otherymp]
    @test !(Γ0_histy[eq_mp, eq_other]  ≈ Γ0_alt_smooth_alt[eq_mp, eq_other])
    @test !(Γ1_histy[eq_mp, eq_other]  ≈ Γ1_alt_smooth_alt[eq_mp, eq_other])

    @test Γ0_histy[eq_othery, eq_otherymp] ≈ Γ0_alt_flex_ait[eq_othery, eq_otherymp]
    @test Γ1_histy[eq_othery, eq_otherymp] ≈ Γ1_alt_flex_ait[eq_othery, eq_otherymp]
    @test !(Γ0_histy[eq_mp, eq_other]  ≈ Γ0_alt_flex_ait[eq_mp, eq_other])
    @test !(Γ1_histy[eq_mp, eq_other]  ≈ Γ1_alt_flex_ait[eq_mp, eq_other])

    @test Γ0_histy[eq_othery, eq_otherymp] ≈ Γ0_alt_ait[eq_othery, eq_otherymp]
    @test Γ1_histy[eq_othery, eq_otherymp] ≈ Γ1_alt_ait[eq_othery, eq_otherymp]
    @test !(Γ0_histy[eq_mp, eq_other]  ≈ Γ0_alt_ait[eq_mp, eq_othery])
    @test !(Γ1_histy[eq_mp, eq_other]  ≈ Γ1_alt_ait[eq_mp, eq_othery])

    @test Γ0_histy[eq_othery, eq_otherymp] ≈ Γ0_alt_ngdp[eq_othery, eq_otherymp]
    @test Γ1_histy[eq_othery, eq_otherymp] ≈ Γ1_alt_ngdp[eq_othery, eq_otherymp]
    @test !(Γ0_histy[eq_mp, eq_other]  ≈ Γ0_alt_ngdp[eq_mp, eq_othery])
    @test !(Γ1_histy[eq_mp, eq_other]  ≈ Γ1_alt_ngdp[eq_mp, eq_othery])

    @test Γ0_histy[eq_othery, eq_otherymp] ≈ Γ0_alt_rw[eq_othery, eq_otherymp]
    @test Γ1_histy[eq_othery, eq_otherymp] ≈ Γ1_alt_rw[eq_othery, eq_otherymp]
    @test !(Γ0_histy[eq_mp, eq_other]  ≈ Γ0_alt_rw[eq_mp, eq_othery])
    @test !(Γ1_histy[eq_mp, eq_other]  ≈ Γ1_alt_rw[eq_mp, eq_othery])

    @test Γ0_histy[eq_othery, eq_otherymp] ≈ Γ0_alt_rw_zero[eq_othery, eq_otherymp]
    @test Γ1_histy[eq_othery, eq_otherymp] ≈ Γ1_alt_rw_zero[eq_othery, eq_otherymp]
    @test !(Γ0_histy[eq_mp, eq_other]  ≈ Γ0_alt_rw_zero[eq_mp, eq_othery])
    @test !(Γ1_histy[eq_mp, eq_other]  ≈ Γ1_alt_rw_zero[eq_mp, eq_othery])
end

# Checking if setting variables works
## These values are random and not actually good choices per se
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
