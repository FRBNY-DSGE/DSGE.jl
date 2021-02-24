using DSGE, Test, JLD2, ModelConstructors, Dates
path = dirname(@__FILE__)

# Set up arguments
m = AnSchorfheide(testing = true)
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m <= Setting(:forecast_horizons, 1)

system = JLD2.jldopen("$path/../reference/forecast_args.jld2","r") do file
    read(file, "system")
end
z0 = zeros(n_states_augmented(m))

# Read expected output
exp_states, exp_obs, exp_pseudo, exp_shocks =
    JLD2.jldopen("$path/../reference/forecast_out.jld2", "r") do file
        read(file, "exp_states"),
        read(file, "exp_obs"),
        read(file, "exp_pseudo"),
        read(file, "exp_shocks")
    end

# Without shocks
global states, obs, pseudo, shocks = forecast(m, system, z0; draw_shocks = false)

@testset "Testing forecasting without drawing shocks" begin
    @test @test_matrix_approx_eq exp_states states
    @test @test_matrix_approx_eq exp_obs    obs
    @test @test_matrix_approx_eq exp_pseudo pseudo
    @test @test_matrix_approx_eq exp_shocks shocks
end

# Supplying shocks
global states, obs, pseudo, shocks = forecast(m, system, z0; shocks = shocks)

@testset "Testing forecasting with pre-supplied shocks" begin
    @test @test_matrix_approx_eq exp_states states
    @test @test_matrix_approx_eq exp_obs    obs
    @test @test_matrix_approx_eq exp_pseudo pseudo
    @test @test_matrix_approx_eq exp_shocks shocks
end

# Draw normally distributed shocks
global states, obs, pseudo, shocks = forecast(m, system, z0; draw_shocks = true)

# Draw t-distributed shocks
m <= Setting(:forecast_tdist_shocks, true)
global states, obs, pseudo, shocks = forecast(m, system, z0; draw_shocks = true)
m <= Setting(:forecast_tdist_shocks, false)

# Enforce ZLB
ind_r = m.observables[:obs_nominalrate]
ind_r_sh = m.exogenous_shocks[:rm_sh]
zlb_value = forecast_zlb_value(m)
shocks = zeros(n_shocks_exogenous(m), forecast_horizons(m))
shocks[ind_r_sh, :] .= -10.

@testset "Ensure valid forecasting at the ZLB" begin
    global states, obs, pseudo, shocks = forecast(m, system, z0; shocks = shocks)
    @test all(x -> x < zlb_value, obs[ind_r, :])

    global states, obs, pseudo, shocks = forecast(m, system, z0; shocks = shocks, enforce_zlb = true)
    @test all(x -> abs(x - zlb_value) < 0.01, obs[ind_r, :])
    @test all(x -> x != -10.,                 shocks[ind_r_sh, :])
end

@testset "Test that forecasting and populating shocks under alternative policy works" begin
    m <= Setting(:regime_switching, true)
    setup_permanent_altpol!(m, AltPolicy(:taylor93, eqcond, solve))
    @test typeof(forecast(m, system, z0; draw_shocks = true)) == NTuple{4, Array{Float64, 2}}
    delete!(DSGE.get_settings(m), :regime_switching)
end

@testset "Enforce ZLB as a temporary alternative policy" begin

    # Set up model for forecast and permanent NGDP
    m = Model1002("ss10"; custom_settings = Dict{Symbol, Setting}(:add_altpolicy_pgap => Setting(:add_altpolicy_pgap, true)))
    m <= Setting(:date_forecast_start,  Date(2020, 6, 30))
    m <= Setting(:date_conditional_end, Date(2020, 6, 30))
    m <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m),
                                                2 => Date(2020, 3, 31),
                                                3 => Date(2020, 6, 30)))
    m <= Setting(:forecast_horizons, 12)
    shocks = zeros(n_shocks_exogenous(m), forecast_horizons(m))
    shocks[m.exogenous_shocks[:b_sh], 1] = -0.6 # suppose massive negative shock to spreads -> MP should drop!
    m <= Setting(:gensys2, false)
    m <= Setting(:regime_switching, true)
    setup_regime_switching_inds!(m)
    m <= Setting(:pgap_value, 12.)
    m <= Setting(:pgap_type, :ngdp)
    setup_permanent_altpol!(m, DSGE.ngdp(); cond_type = :none)
    system = compute_system(m)
    z0 = zeros(n_states_augmented(m))

    # First test with unconditional
    ngdp_states, ngdp_obs, ngdp_pseudo, _ = forecast(m, system, z0; shocks = shocks, cond_type = :none)
    @test !all(ngdp_obs[m.observables[:obs_nominalrate], :] .> -1e-14)
    ngdp_states, ngdp_obs, ngdp_pseudo = forecast(m, z0, ngdp_states, ngdp_obs, ngdp_pseudo, shocks; cond_type = :none)
    @test all(ngdp_obs[m.observables[:obs_nominalrate], :] .> -1e-14)

    # Now test a conditional forecast with regime switching in the forecast
    m <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m),
                                                2 => Date(2020, 3, 31),
                                                3 => Date(2020, 6, 30),
                                                4 => Date(2020, 9, 30)))
    m <= Setting(:forecast_horizons, 12)
    shocks = zeros(n_shocks_exogenous(m), forecast_horizons(m; cond_type = :full))
    shocks[m.exogenous_shocks[:b_sh], 1] = -0.6 # suppose massive negative shock to spreads -> MP should drop!
    m <= Setting(:replace_eqcond, false)
    m <= Setting(:gensys2, false)
    m <= Setting(:regime_switching, true)
    setup_regime_switching_inds!(m)
    m <= Setting(:pgap_value, 12.)
    m <= Setting(:pgap_type, :ngdp)
    setup_permanent_altpol!(m, DSGE.ngdp(); cond_type = :full)
    system = compute_system(m)
    z0 = zeros(n_states_augmented(m))

    ngdp_states, ngdp_obs, ngdp_pseudo, _ = forecast(m, system, z0; shocks = shocks, cond_type = :full)
    @test !all(ngdp_obs[m.observables[:obs_nominalrate], :] .> -1e-14)
    ngdp_states, ngdp_obs, ngdp_pseudo = forecast(m, z0, ngdp_states, ngdp_obs, ngdp_pseudo, shocks; cond_type = :full)
    @test all(ngdp_obs[m.observables[:obs_nominalrate], :] .> -1e-14)
end

nothing
