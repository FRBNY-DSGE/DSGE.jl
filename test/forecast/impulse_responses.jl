path = dirname(@__FILE__)

# Set up arguments
global m = AnSchorfheide(testing = true)

global system = JLD2.jldopen("$path/../reference/forecast_args.jld2","r") do file
    read(file, "system")
end

# Run impulse responses
states, obs, pseudo = impulse_responses(m, system)

# Compare to expected output
exp_states, exp_obs, exp_pseudo =
    JLD2.jldopen("$path/../reference/impulse_responses_out.jld2", "r") do file
        read(file, "exp_states"), read(file, "exp_obs"), read(file, "exp_pseudo")
    end

@testset "Compare irfs to expected output" begin
    @test @test_matrix_approx_eq exp_states states
    @test @test_matrix_approx_eq exp_obs    obs
    @test @test_matrix_approx_eq exp_pseudo pseudo
end

# Test impulse response methods for:
# (1) Choosing a set of shocks and their shock values
# (2) Choosing a shock, a state, and the value the state should hit in period 1 based on
# an impulse from that shock at time 1.
# (3) Same as (2) but with an observable.

# Setup
endo   = m.endogenous_states
obsvs  = m.observables
exo    = m.exogenous_shocks
horizon = 10

# Testing both for exact output match and for certain values to match to be
# explicit about the expected outcome of the impulse response

# Shock
shock_names  = [:z_sh, :g_sh]
shock_values = [1., 1.]
states, obs, pseudo = impulse_responses(m, system, horizon, shock_names, shock_values)

# Explicit checks
@testset "Compare shock-set irfs to the states they affect 1-for-1" begin
    @test states[endo[:z_t], 1, exo[:z_sh]] ≈ 1.
    @test states[endo[:y_t], 1, exo[:g_sh]] ≈ 1.
    @test states[endo[:g_t], 1, exo[:g_sh]] ≈ 1.
end

exp_states_shockset, exp_obs_shockset, exp_pseudo_shockset =
    jldopen("$path/../reference/impulse_responses_out.jld2", "r") do file
        read(file, "exp_states_shockset"), read(file, "exp_obs_shockset"), read(file, "exp_pseudo_shockset")
    end

@testset "Checking shockset irfs" begin
    @test states ≈ exp_states_shockset
    @test obs    ≈ exp_obs_shockset
    @test pseudo ≈ exp_pseudo_shockset
end

# State
shock_name = :z_sh
state_name   = :y_t
state_value  = 1.
states, obs, pseudo = impulse_responses(m, system, horizon, shock_name, state_name, state_value)

# Explicit checks
@testset "Check state irf initial period is exactly the specified state value" begin
    @test states[endo[:y_t], 1, exo[:z_sh]] ≈ 1.
end

exp_states_shockstates, exp_obs_shockstates, exp_pseudo_shockstates =
    jldopen("$path/../reference/impulse_responses_out.jld2", "r") do file
        read(file, "exp_states_shockstate"), read(file, "exp_obs_shockstate"), read(file, "exp_pseudo_shockstate")
    end

@testset "Checking shockstate irfs" begin
    @test states ≈ exp_states_shockstates
    @test obs    ≈ exp_obs_shockstates
    @test pseudo ≈ exp_pseudo_shockstates
end

# Observable
shock_name = :z_sh
obs_name   = :obs_gdp
obs_value  = 1.
states, obs, pseudo = impulse_responses(m, system, horizon, shock_name, obs_name, obs_value)

# Explicit checks
@testset "Check obs irf initial period is exactly the specified obs value" begin
    @test obs[obsvs[:obs_gdp], 1, exo[:z_sh]] ≈ 1.
end

exp_states_shockobs, exp_obs_shockobs, exp_pseudo_shockobs =
    jldopen("$path/../reference/impulse_responses_out.jld2", "r") do file
        read(file, "exp_states_shockobs"), read(file, "exp_obs_shockobs"), read(file, "exp_pseudo_shockobs")
    end

@testset "Checking shockobs irfs" begin
    @test states ≈ exp_states_shockobs
    @test obs    ≈ exp_obs_shockobs
    @test pseudo ≈ exp_pseudo_shockobs
end

# Test impulse response methods for:
# (1) Choosing a set of shocks and their shock values
# (2) Choosing a shock, a state, and the value the state should hit in period 1 based on
# an impulse from that shock at time 1.
# (3) Same as (2) but with an observable.

# Setup
endo   = m.endogenous_states
obsvs  = m.observables
exo    = m.exogenous_shocks
horizon = 10

# Testing both for exact output match and for certain values to match to be
# explicit about the expected outcome of the impulse response

# Shock
shock_names  = [:z_sh, :g_sh]
shock_values = [1., 1.]
states, obs, pseudo = impulse_responses(m, system, horizon, shock_names, shock_values)

# Explicit checks
@testset "Compare shock-set irfs to the states they affect 1-for-1" begin
    @test states[endo[:z_t], 1, exo[:z_sh]] ≈ 1.
    @test states[endo[:y_t], 1, exo[:g_sh]] ≈ 1.
    @test states[endo[:g_t], 1, exo[:g_sh]] ≈ 1.
end

exp_states_shockset, exp_obs_shockset, exp_pseudo_shockset =
    jldopen("$path/../reference/impulse_responses_out.jld2", "r") do file
        read(file, "exp_states_shockset"), read(file, "exp_obs_shockset"), read(file, "exp_pseudo_shockset")
    end

@testset "Checking shockset irfs" begin
    @test states ≈ exp_states_shockset
    @test obs    ≈ exp_obs_shockset
    @test pseudo ≈ exp_pseudo_shockset
end

# State
shock_name = :z_sh
state_name   = :y_t
state_value  = 1.
states, obs, pseudo = impulse_responses(m, system, horizon, shock_name, state_name, state_value)

# Explicit checks
@testset "Check state irf initial period is exactly the specified state value" begin
    @test states[endo[:y_t], 1, exo[:z_sh]] ≈ 1.
end

exp_states_shockstates, exp_obs_shockstates, exp_pseudo_shockstates =
    jldopen("$path/../reference/impulse_responses_out.jld2", "r") do file
        read(file, "exp_states_shockstate"), read(file, "exp_obs_shockstate"), read(file, "exp_pseudo_shockstate")
    end

@testset "Checking shockstate irfs" begin
    @test states ≈ exp_states_shockstates
    @test obs    ≈ exp_obs_shockstates
    @test pseudo ≈ exp_pseudo_shockstates
end

# Observable
shock_name = :z_sh
obs_name   = :obs_gdp
obs_value  = 1.
states, obs, pseudo = impulse_responses(m, system, horizon, shock_name, obs_name, obs_value)

# Explicit checks
@testset "Check obs irf initial period is exactly the specified obs value" begin
    @test obs[obsvs[:obs_gdp], 1, exo[:z_sh]] ≈ 1.
end

exp_states_shockobs, exp_obs_shockobs, exp_pseudo_shockobs =
    jldopen("$path/../reference/impulse_responses_out.jld2", "r") do file
        read(file, "exp_states_shockobs"), read(file, "exp_obs_shockobs"), read(file, "exp_pseudo_shockobs")
    end

@testset "Checking shockobs irfs" begin
    @test states ≈ exp_states_shockobs
    @test obs    ≈ exp_obs_shockobs
    @test pseudo ≈ exp_pseudo_shockobs
end

# Test impulse response method for computing
# a short-run Cholesky-identified shock
obs_shock = zeros(n_observables(m))
obs_shock[1] = 1.
states_chol, obs_chol, pseudo_chol, struct_shock =
    impulse_responses(system, horizon, Matrix{Float64}(I, n_observables(m), n_observables(m)),
                      obs_shock, flip_shocks = false, get_shocks = true)
states_chol1, obs_chol1, pseudo_chol1 =
    impulse_responses(system, horizon, Matrix{Float64}(I, n_observables(m), n_observables(m)),
                      obs_shock, flip_shocks = false, get_shocks = false)
states_chol2, obs_chol2, pseudo_chol2 =
    impulse_responses(system, horizon, Matrix{Float64}(I, n_observables(m), n_observables(m)),
                      flip_shocks = false, get_shocks = false)
states_chol3, obs_chol3, pseudo_chol3 =
    impulse_responses(system, horizon, Matrix{Float64}(I, n_observables(m), n_observables(m)),
                      flip_shocks = true, get_shocks = false)

exp_states_chol, exp_obs_chol, exp_pseudo_chol, exp_struct_shock =
    JLD2.jldopen("$path/../reference/impulse_responses_out.jld2", "r") do file
        read(file, "exp_states_chol"), read(file, "exp_obs_chol"), read(file, "exp_pseudo_chol"),
        read(file, "exp_struct_shock")
    end

@testset "Compare irfs to expected output for a short-run Cholesky-identified shock" begin
    @test @test_matrix_approx_eq states_chol1 states_chol
    @test @test_matrix_approx_eq obs_chol1 obs_chol
    @test @test_matrix_approx_eq pseudo_chol1  pseudo_chol
    @test @test_matrix_approx_eq states_chol2 states_chol
    @test @test_matrix_approx_eq obs_chol2 obs_chol
    @test @test_matrix_approx_eq pseudo_chol2  pseudo_chol
    @test @test_matrix_approx_eq states_chol3 -states_chol
    @test @test_matrix_approx_eq obs_chol3 -obs_chol
    @test @test_matrix_approx_eq pseudo_chol3  -pseudo_chol
    @test @test_matrix_approx_eq states_chol  exp_states_chol
    @test @test_matrix_approx_eq obs_chol  exp_obs_chol
    @test @test_matrix_approx_eq pseudo_chol  exp_pseudo_chol
    @test @test_matrix_approx_eq struct_shock exp_struct_shock
    @test @test_matrix_approx_eq -system[:RRR]*struct_shock states_chol[:, 1]
    @test @test_matrix_approx_eq -system[:ZZ]*system[:RRR]*struct_shock obs_chol[:, 1]
end


# Test impulse response method for computing
# a long-run Cholesky-identified shock
obs_shock = zeros(n_observables(m))
obs_shock[1] = 1.
states_chol, obs_chol, pseudo_chol, struct_shock =
    impulse_responses(system, horizon, Matrix{Float64}(I, n_observables(m), n_observables(m)),
                      obs_shock, flip_shocks = false, get_shocks = true,
                      restriction = :long_run)
states_chol1, obs_chol1, pseudo_chol1 =
    impulse_responses(system, horizon, Matrix{Float64}(I, n_observables(m), n_observables(m)),
                      obs_shock, flip_shocks = false, get_shocks = false,
                      restriction = :long_run)
states_chol2, obs_chol2, pseudo_chol2 =
    impulse_responses(system, horizon, Matrix{Float64}(I, n_observables(m), n_observables(m)),
                      flip_shocks = false, get_shocks = false,
                      restriction = :long_run)
states_chol3, obs_chol3, pseudo_chol3 =
    impulse_responses(system, horizon, Matrix{Float64}(I, n_observables(m), n_observables(m)),
                      flip_shocks = true, get_shocks = false, restriction = :long_run)

exp_states_chol, exp_obs_chol, exp_pseudo_chol, exp_struct_shock =
    JLD2.jldopen("$path/../reference/impulse_responses_out.jld2", "r") do file
        read(file, "exp_states_chollr"), read(file, "exp_obs_chollr"), read(file, "exp_pseudo_chollr"),
        read(file, "exp_struct_shocklr")
    end

@testset "Compare irfs to expected output for a long-run Cholesky-identified shock" begin
    @test @test_matrix_approx_eq states_chol1 states_chol
    @test @test_matrix_approx_eq obs_chol1 obs_chol
    @test @test_matrix_approx_eq pseudo_chol1  pseudo_chol
    @test @test_matrix_approx_eq states_chol2 states_chol
    @test @test_matrix_approx_eq obs_chol2 obs_chol
    @test @test_matrix_approx_eq pseudo_chol2  pseudo_chol
    @test @test_matrix_approx_eq states_chol3 -states_chol
    @test @test_matrix_approx_eq obs_chol3 -obs_chol
    @test @test_matrix_approx_eq pseudo_chol3  -pseudo_chol
    @test @test_matrix_approx_eq states_chol  exp_states_chol
    @test @test_matrix_approx_eq obs_chol  exp_obs_chol
    @test @test_matrix_approx_eq pseudo_chol  exp_pseudo_chol
    @test @test_matrix_approx_eq struct_shock exp_struct_shock
    @test @test_matrix_approx_eq -system[:RRR]*struct_shock states_chol[:, 1]
    @test @test_matrix_approx_eq -system[:ZZ]*system[:RRR]*struct_shock obs_chol[:, 1]
end

# Test impulse response method for computing
# a shock to maximizes business-cycle variance
states_bc1, obs_bc1, pseudo_bc1 = impulse_responses(system, 4, (2 * π / 32, 2 * π / 6), 2)
states_bc2, _ = impulse_responses(system, 4, (2 * π / 32, 2 * π / 6), 2,
                                                    flip_shocks = true)
exp_states_bc, exp_obs_bc, exp_pseudo_bc =
    JLD2.jldopen("$path/../reference/impulse_responses_out.jld2", "r") do file
        read(file, "exp_states_bc"), read(file, "exp_obs_bc"), read(file, "exp_pseudo_bc")
    end

@testset "Compare irfs to expected output for shock maximizing business cycle variance of an observable" begin
    @test @test_matrix_approx_eq states_bc1 -states_bc2
    @test @test_matrix_approx_eq states_bc1  exp_states_bc
    @test @test_matrix_approx_eq obs_bc1  exp_obs_bc
    @test @test_matrix_approx_eq pseudo_bc1  exp_pseudo_bc
end


@testset "Check impulse responses to pegged interest rate" begin
    global m = SmetsWouters(
        custom_settings = Dict{Symbol, Setting}(:n_anticipated_shocks => Setting(:n_anticipated_shocks, 8)))
    global system = compute_system(m)
    for hor=1:8
        DSGE.impulse_responses_peg(m, system, horizon, H = hor, peg = :all_periods, real_rate = false)
        DSGE.impulse_responses_peg(m, system, horizon, H = hor, peg = :all_periods, real_rate = true)
        DSGE.impulse_responses_peg(m, system, horizon, H = hor, peg = :some_periods, real_rate = false)
        DSGE.impulse_responses_peg(m, system, horizon, H = hor, peg = :some_periods, real_rate = true)
    end
end


nothing
