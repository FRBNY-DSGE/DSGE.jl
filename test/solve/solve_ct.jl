using DSGE, JLD2, DataStructures
import Test: @test, @testset
import DSGE: @test_matrix_approx_eq

# Load true values
true_inverse_basis = load("../../test_outputs/solve/solve.jld2", "inverse_basis")
G1 = load("../../test_outputs/solve/solve.jld2", "G1")
impact = load("../../test_outputs/solve/solve.jld2", "impact")

# Test
m = KrusellSmith()
### Steady state

# Load test values
params = load("../../test_outputs/krusell_smith/saved_params.jld2", "params")
grids = load("../../test_outputs/krusell_smith/saved_params.jld2", "grids")
init_params = load("../../test_outputs/krusell_smith/saved_params.jld2", "init_params")
approx_params = load("../../test_outputs/krusell_smith/saved_params.jld2", "approx_params")
vars_ss_mat = load("../../test_outputs/krusell_smith/vars_ss.jld2")
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
update!(m.settings[:I], Setting(:I, 100))
update!(m.settings[:amin], Setting(:amin, 0.))
update!(m.settings[:amax], Setting(:amax, 100.))
update!(m.settings[:a], Setting(:a, collect(linspace(0., 100., 100))))
update!(m.settings[:da], Setting(:da, 100./(100 - 1)))
update!(m.settings[:aaa], Setting(:aaa, reshape(get_setting(m, :aa), 2 * 100, 1)))
update!(m.settings[:J], Setting(:J, 2))
update!(m.settings[:z], Setting(:z, [0., 1.]))
update!(m.settings[:dz], Setting(:dz, 1.))
update!(m.settings[:zz], Setting(:zz, ones(100, 1) * get_setting(m, :z)'))
update!(m.settings[:zzz], Setting(:zzz, reshape(get_setting(m, :zz), 2 * 100, 1)))
update!(m.settings[:n_jump_vars], Setting(:n_jump_vars, length(m.endogenous_states[:value_function])))
update!(m.settings[:n_state_vars], Setting(:n_state_vars, length(m.endogenous_states[:distribution]) + n_shocks_exogenous(m)))
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
update!(m.settings[:knots_dict], Setting(:knots_dict, Dict(1 => collect(linspace(get_setting(m, :amin), get_setting(m, :amax), 12 - 1)))))
update!(m.settings[:spline_grid], Setting(:spline_grid, get_setting(m, :a)))
update!(m.settings[:n_prior], Setting(:n_prior, 1))
update!(m.settings[:n_post], Setting(:n_post, 2))

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
