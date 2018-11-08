using JLD
import DataStructures: OrderedDict
import Test: @test, @testset
import DSGE: @test_matrix_approx_eq
using DSGE

### Model
m = KrusellSmithCT()

### Model indices

# Endogenous states
endo = m.endogenous_states
@testset "Model indices" begin
    @test n_states(m) == length(endo[:value_function]) + length(endo[:distribution]) + length(endo[:log_aggregate_tfp]) + 6
    @test get_setting(m, :n_jump_vars) == length(endo[:value_function])
    @test get_setting(m, :n_state_vars) == length(endo[:distribution]) + length(endo[:log_aggregate_tfp]) - get_setting(m, :n_state_vars_unreduce)

    # Exogenous shocks
    @test n_shocks_exogenous(m) == 1

    # Expectation shocks
    @test n_shocks_expectational(m) == length(endo[:value_function])
end

### Steady state
#=
# Load test values
# TODO: These file paths are out of date
params = load("../../../test_outputs/krusell_smith/saved_params.jld", "params")
grids = load("../../../test_outputs/krusell_smith/saved_params.jld", "grids")
init_params = load("../../../test_outputs/krusell_smith/saved_params.jld", "init_params")
approx_params = load("../../../test_outputs/krusell_smith/saved_params.jld", "approx_params")
vars_ss_mat = load("../../../test_outputs/krusell_smith/vars_ss.jld")
vars_ss_mat = vars_ss_mat["varsSS"]

# Update parameters
m[:γ] = params[:γ]
m[:ρ] = params[:ρ]
m[:δ] = params[:δ]
m[:α] = params[:α]
m[:σ_tfp] = params[:σ_tfp]
m[:ρ_tfp] = params[:ρ_tfp]
m[:λ1] = params[:λ1]
m[:λ2] = params[:λ2]
m[:μ] = params[:μ]
m[:τ] = params[:τ]
=#
# Update settings
update!(m.settings[:I], Setting(:I, 100))
update!(m.settings[:amin], Setting(:amin, 0.))
update!(m.settings[:amax], Setting(:amax, 100.))
update!(m.settings[:a], Setting(:a, collect(range(0., stop=100., length=100))))
update!(m.settings[:da], Setting(:da, 100 ./ (100 - 1)))
update!(m.settings[:aaa], Setting(:aaa, reshape(get_setting(m, :aa), 2 * 100, 1)))
update!(m.settings[:J], Setting(:J, 2))
update!(m.settings[:z], Setting(:z, [0.; 1.]))
update!(m.settings[:dz], Setting(:dz, 1.))
update!(m.settings[:zz], Setting(:zz, ones(100, 1) * get_setting(m, :z)'))
update!(m.settings[:zzz], Setting(:zzz, reshape(get_setting(m, :zz), 2 * 100, 1)))
endo = m.endogenous_states
update!(m.settings[:n_jump_vars], Setting(:n_jump_vars, length(endo[:value_function])))
update!(m.settings[:n_state_vars], Setting(:n_state_vars, length(endo[:distribution]) .+ n_shocks_exogenous(m)))
update!(m.settings[:n_state_vars_unreduce], Setting(:n_state_vars_unreduce, 0))
update!(m.settings[:rmax], Setting(:rmin, .0001))
update!(m.settings[:rmax], Setting(:rmax, m[:ρ].value))
update!(m.settings[:r0], Setting(:r0, .005))
update!(m.settings[:maxit_HJB], Setting(:maxit_HJB, 100))
update!(m.settings[:crit_HJB], Setting(:crit_HJB, 1e-6))
update!(m.settings[:Δ_HJB], Setting(:Δ_HJB, 1e4))
update!(m.settings[:Ir], Setting(:Ir, 100))
update!(m.settings[:crit_S], Setting(:crit_S, 1e-5))
update!(m.settings[:reduce_state_vars], Setting(:reduce_state_vars, true))
update!(m.settings[:reduce_v], Setting(:reduce_v, true))
update!(m.settings[:krylov_dim], Setting(:krylov_dim, 5))
update!(m.settings[:n_knots], Setting(:n_knots, 12))
update!(m.settings[:c_power], Setting(:c_power, 7))
update!(m.settings[:knots_dict], Setting(:knots_dict, Dict(1 => collect(range(get_setting(m, :amin),
                                                                              stop=get_setting(m, :amax),
                                                                              length=(12 - 1))))))
update!(m.settings[:spline_grid], Setting(:spline_grid, get_setting(m, :a)))
update!(m.settings[:n_prior], Setting(:n_prior, 1))
update!(m.settings[:n_post], Setting(:n_post, 2))
#=
# compare steady state values
steadystate!(m)
@testset "Steady State" begin
    @test @test_matrix_approx_eq vec(m[:V_ss].value) vec(vars_ss_mat[:VSS])
    @test @test_matrix_approx_eq vec(m[:gg_ss].value) vec(vars_ss_mat[:ggSS])
    @test m[:K_ss].value ≈ vars_ss_mat[:KSS] # allow for machine error
    @test m[:r_ss].value ≈ vars_ss_mat[:rSS]
    @test m[:w_ss].value ≈ vars_ss_mat[:wSS]
    @test m[:Y_ss].value ≈ vars_ss_mat[:YSS]
    @test m[:C_ss].value ≈ vars_ss_mat[:CSS]
    @test m[:I_ss].value ≈ vars_ss_mat[:ISS]
end
=#
### Equilibrium conditions
Γ0, Γ1, Ψ, Π, C = eqcond(m)

# Transition and measurement equations
TTT, RRR, CCC, ~ = solve(m)
# measurement equation here

# Matrices have expected dimensions
@testset "Check that outputted matrices from solve have the right dimensions" begin
    n_vars = n_states(m)
    @test size(Γ0) == (n_vars, n_vars)
    @test size(Γ1) == (n_vars, n_vars)
    @test size(vec(C)) == (n_vars,)
    @test size(Ψ) == (n_vars, 1)
    @test size(Π) == (n_vars, n_shocks_expectational(m))
    if get_setting(m, :reduce_v)
        @test size(TTT) == (get_setting(m, :n_splined) + get_setting(m, :n_state_vars_red), get_setting(m, :n_splined) + get_setting(m, :n_state_vars_red))
        @test size(RRR) == (get_setting(m, :n_splined) + get_setting(m, :n_state_vars_red), 1)
        @test size(vec(CCC)) == (get_setting(m, :n_splined) + get_setting(m, :n_state_vars_red),)
    elseif get_setting(m, :reduce_state_vars)
        @test size(TTT) == (get_setting(m, :n_jump_vars) + get_setting(m, :n_state_vars_red), get_setting(m, :n_state_vars_red))
        @test size(RRR) == (get_setting(m, :n_jump_vars) + get_setting(m, :n_state_vars_red), 1)
        @test size(vec(CCC)) == (get_setting(m, :n_jump_vars) + get_setting(m, :n_state_vars_red),)
    else
        @test size(TTT) == (n_vars, n_vars)
        @test size(RRR) == (n_vars, 1)
        @test size(vec(CCC)) == (n_vars,)
    end
end
