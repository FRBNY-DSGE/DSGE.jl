using DSGE, JLD, DataStructures
import Test: @test, @testset
import DSGE: @test_matrix_approx_eq

# Load true values
true_inverse_basis = load("../../test_outputs/solve/solve.jld", "inverse_basis")
G1 = load("../../test_outputs/solve/solve.jld", "G1")
impact = load("../../test_outputs/solve/solve.jld", "impact")

# Test
m = KrusellSmith()
### Steady state

# Load test values
params = load("../../test_outputs/krusell_smith/saved_params.jld", "params")
grids = load("../../test_outputs/krusell_smith/saved_params.jld", "grids")
init_params = load("../../test_outputs/krusell_smith/saved_params.jld", "init_params")
approx_params = load("../../test_outputs/krusell_smith/saved_params.jld", "approx_params")
vars_ss_mat = load("../../test_outputs/krusell_smith/vars_ss.jld")
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

# Update settings
EHANK.update!(m.settings[:I], Setting(:I, 100))
EHANK.update!(m.settings[:amin], Setting(:amin, 0.))
EHANK.update!(m.settings[:amax], Setting(:amax, 100.))
EHANK.update!(m.settings[:a], Setting(:a, collect(linspace(0., 100., 100))))
EHANK.update!(m.settings[:da], Setting(:da, 100./(100 - 1)))
EHANK.update!(m.settings[:aaa], Setting(:aaa, reshape(get_setting(m, :aa), 2 * 100, 1)))
EHANK.update!(m.settings[:J], Setting(:J, 2))
EHANK.update!(m.settings[:z], Setting(:z, [0., 1.]))
EHANK.update!(m.settings[:dz], Setting(:dz, 1.))
EHANK.update!(m.settings[:zz], Setting(:zz, ones(100, 1) * get_setting(m, :z)'))
EHANK.update!(m.settings[:zzz], Setting(:zzz, reshape(get_setting(m, :zz), 2 * 100, 1)))
EHANK.update!(m.settings[:n_jump_vars], Setting(:n_jump_vars, length(m.endogenous_states[:value_function])))
EHANK.update!(m.settings[:n_state_vars], Setting(:n_state_vars, length(m.endogenous_states[:distribution]) + n_shocks_exogenous(m)))
EHANK.update!(m.settings[:n_state_vars_unreduce], Setting(:n_state_vars_unreduce, 0))
EHANK.update!(m.settings[:rmax], Setting(:rmin, .0001))
EHANK.update!(m.settings[:rmax], Setting(:rmax, m[:ρ].value))
EHANK.update!(m.settings[:r0], Setting(:r0, .005))
EHANK.update!(m.settings[:maxit_HJB], Setting(:maxit_HJB, 100))
EHANK.update!(m.settings[:crit_HJB], Setting(:crit_HJB, 1e-6))
EHANK.update!(m.settings[:Δ_HJB], Setting(:Δ_HJB, 1e4))
EHANK.update!(m.settings[:Ir], Setting(:Ir, 100))
EHANK.update!(m.settings[:crit_S], Setting(:crit_S, 1e-5))
EHANK.update!(m.settings[:reduce_state_vars], Setting(:reduce_state_vars, true))
EHANK.update!(m.settings[:reduce_v], Setting(:reduce_v, true))
EHANK.update!(m.settings[:krylov_dim], Setting(:krylov_dim, 5))
EHANK.update!(m.settings[:n_knots], Setting(:n_knots, 12))
EHANK.update!(m.settings[:c_power], Setting(:c_power, 7))
EHANK.update!(m.settings[:knots_dict], Setting(:knots_dict, Dict(1 => collect(linspace(get_setting(m, :amin), get_setting(m, :amax), 12 - 1)))))
EHANK.update!(m.settings[:spline_grid], Setting(:spline_grid, get_setting(m, :a)))
EHANK.update!(m.settings[:n_prior], Setting(:n_prior, 1))
EHANK.update!(m.settings[:n_post], Setting(:n_post, 2))

# Perform solve
steadystate!(m)
TTT, RRR, ~, inverse_basis, ~ = solve(m; sparse_mat = false)
@testset "Simulating Impulse Responses (non-sparse)" begin
    @test @test_matrix_approx_eq full(TTT) full(G1)
    @test @test_matrix_approx_eq full(RRR) full(impact)
    @test @test_matrix_approx_eq full(true_inverse_basis) full(inverse_basis)
end

steadystate!(m)
TTT, RRR, ~, inverse_basis, ~ = solve(m; sparse_mat = true)
@testset "Simulating Impulse Responses (sparse)" begin
    @test @test_matrix_approx_eq full(TTT) full(G1)
    @test @test_matrix_approx_eq full(RRR) full(impact)
    @test @test_matrix_approx_eq full(true_inverse_basis) full(inverse_basis)
end
