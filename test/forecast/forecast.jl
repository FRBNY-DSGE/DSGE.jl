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

nothing
