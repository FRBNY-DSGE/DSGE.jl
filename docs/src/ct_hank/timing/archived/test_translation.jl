using JLD, DataStructures, MAT
using ForwardDiff, Plots
using Compat, BenchmarkTools
gr()

### The test file for comparing Julia and Matlab HANK codes

include("set_parameters.jl")
include("util.jl")
include("helpers.jl")
include("compute_steady_state.jl")
include("equilibrium_conditions.jl")
include("solve.jl")
include("simulate.jl")
include("simulate_irfs.jl")
include("gensysct.jl")
include("krylov_reduction.jl")

ReduceDistribution = 1
reduceV = 1
ReduceDist_hors = 20
DisplayLev = 1
check_consistency = 1
atol = 1e-10
krylov_dim = 20
run_timing_tests = true

params, approx_params, grids, I, J = set_parameters()
println("Parameters set.")

n_v = Int64(I*J + 1) # The number of jump variables (value function + inflation)
n_g = Int64(I*J)     # The number of endogenous state variables (distribution + monetary policy)
n_p = 5              # The number of static relations: bond-market clearing, labor
                     # market clearing, consumption, output, total assets
n_shocks = 1         # only monetary policy shock is considered
n_exp_errors = n_v
n_vars = n_v + n_g + n_p

## check residual
println("computing steady state...")
vars_SS = compute_steady_state(grids, params, approx_params; DisplayLev = 0)
run_timing_tests && @btime compute_steady_state($grids, $params, $approx_params; DisplayLev = 0);

vars_SS_matlab = matread("saved_outputs/vars_SS.mat")
for var in [:B_SS, :Y_SS, :inflation_SS, :N_SS, :C_SS,
            :V_SS, :T_SS, :profit_SS, :rnom_SS,
            :s_SS,  :w_SS, :c_SS, :h_SS, :g_SS]
    @assert vars_SS_matlab[string(var)] ≈ vars_SS[var] "equality test failed for $var"
end
println("taking derivatives in Julia succeeded")

x = zeros(Float64, 2*n_vars + n_exp_errors + n_shocks)
v_residual = equilibrium_conditions(x, vars_SS, grids, params, n_v, n_g, n_p, niter_hours = 10)

## take derivatives, compare
println("calculating derivatives...")
function f{T<:Real}(x::Vector{T}; vars_ss = deepcopy(vars_SS))
    v_residual = equilibrium_conditions(x, vars_ss, grids, params, n_v, n_g, n_p, niter_hours = 10)
    return vcat(values(v_residual)...)
end
derivs = ForwardDiff.jacobian(f, x)
run_timing_tests && @btime ForwardDiff.jacobian($f, $x);
# derivs_plot = spy(derivs);

derivs_output = matread("saved_outputs/derivs.mat")
derivs_matlab = Matrix(derivs_output["derivs"])

Γ1 = -derivs[:,1:n_vars];
Γ0 = derivs[:,n_vars+1:2*n_vars];
Π = -derivs[:,2*n_vars+1:2*n_vars+n_exp_errors];
Ψ = -derivs[:,2*n_vars+n_exp_errors+1:2*n_vars+n_exp_errors+n_shocks];
C = zeros(n_vars, 1)

pre_reduction_output = matread("saved_outputs/pre_reduction.mat")
g0  = Matrix(pre_reduction_output["g0"]); g1  = Matrix(pre_reduction_output["g1"])
pi  = Matrix(pre_reduction_output["pi"]); psi = Matrix(pre_reduction_output["psi"])

try
    @assert isapprox(g0, Γ0, atol = atol); @assert isapprox(g1, Γ1, atol = atol)
    @assert isapprox(pi, Π, atol = atol); @assert isapprox(psi, Ψ, atol = atol)
    println("taking derivatives in Julia succeeded")
catch
    # For now, just use the read-in Matlab matrices
    println("taking derivatives in Julia failed")
    Γ1 = g1; Γ0 = g0
    Π = pi; Ψ = psi
end


## solve out static constraints and then solve model
println("solving out static constraints...")
Γ0_red, Γ1_red, Ψ_red, Π_red, C_red, redundant_states_inv = solve_static_constraints(Γ0, Γ1, Ψ, Π, C)
run_timing_tests && @btime solve_static_constraints($Γ0, $Γ1, $Ψ, $Π, $C);

reduced_output = matread("saved_outputs/post_reduction.mat")
g0  = Matrix(reduced_output["g0"]); g1  = Matrix(reduced_output["g1"])
pi  = Matrix(reduced_output["pi"]); psi = Matrix(reduced_output["psi"])
try
    @assert isapprox(g0, Γ0_red, atol = atol); @assert isapprox(g1, Γ1_red, atol = atol)
    @assert isapprox(pi, Π_red, atol = atol); @assert isapprox(psi, Ψ_red, atol = atol)
    println("reducing model in Julia succeeded")
catch
    println("reducing model in Julia failed")
    Γ1_red = g1; Γ0_red = g0
    Π_red = pi; Ψ_red = psi
end

println("solving the model...")
T, C, R, ~, ~, ~, ~, eu = gensysct(Γ0_red, Γ1_red, C_red, Ψ_red, Π_red)
run_timing_tests && @btime gensysct($Γ0_red, $Γ1_red, $C_red, $Ψ_red, $Π_red)
gensys_output = matread("saved_outputs/post_gensys.mat")
G1 = Matrix(gensys_output["G1"]); impact = Matrix(gensys_output["impact"])

## simulate IRFs
println("simulating IRFs...")

periods = 100
steps = 400
dt = periods/steps
agg_shock = zeros(1, steps)
agg_shock[1] = 1/sqrt(dt)
transformation = redundant_states_inv

simulated_states = simulate(T, R, periods, steps, agg_shock; method = :implicit)
run_timing_tests && @btime simulate($T, $R, $periods, $steps, $agg_shock; method = :implicit)
simulated_states = (transformation * simulated_states)[cat(1,n_v, n_v+n_g:n_v+n_g+5),:]

if isapprox(T, G1) && isapprox(R, impact)
    println("Gensys in Julia suceeded")
else
    println("Gensys in Julia failed")
end

simulated_states_matlab = simulate(G1, impact, periods, steps, agg_shock; method = :implicit)
simulated_states_matlab = (transformation * simulated_states_matlab)[cat(1,n_v, n_v+n_g:n_v+n_g+5),:]
try
    @assert isapprox(simulated_states, simulated_states_matlab, atol = 1e-3)
    println("IRFs in Julia succeeded (still using Julia's Gensys outputs)")
catch
    println("IRFs in Julia failed (still using Julia's Gensys outputs)")
end

inflation, monetary_shock, wage, consumption, labor_supply, output =
    simulate_irfs(T, R, vars_SS, periods, steps, transformation; DisplayLev = DisplayLev)

## Perform state space and V reduction, then solve model

# Remake state space matrices
Γ1 = -derivs[:,1:n_vars];
Γ0 = derivs[:,n_vars+1:2*n_vars];
Π = -derivs[:,2*n_vars+1:2*n_vars+n_exp_errors];
Ψ = -derivs[:,2*n_vars+n_exp_errors+1:2*n_vars+n_exp_errors+n_shocks];
C = zeros(n_vars, 1)

println("performing Krylov reduction on state space...")

Γ0_kry, Γ1_kry, Ψ_kry, Π_kry, C_kry, basis, basis_inv, ~ = krylov_reduction(Γ0, Γ1, Ψ, Π, C, n_v, n_g, krylov_dim)
run_timing_tests && @btime krylov_reduction($Γ0, $Γ1, $Ψ, $Π, $C, $n_v, $n_g, $krylov_dim);

reduced_output = matread("saved_outputs/krylov_outputs.h5")
g0  = Matrix(reduced_output["g0"]); g1  = Matrix(reduced_output["g1"])
pi  = Matrix(reduced_output["pi"]); psi = Matrix(reduced_output["psi"])

kry_atol = 1e-4
try
    @assert isapprox(g0, Γ0_kry, atol = kry_atol); @assert isapprox(g1, Γ1_kry, atol = kry_atol)
    @assert isapprox(pi, Π_kry, atol = kry_atol); @assert isapprox(psi, Ψ_kry, atol = kry_atol)
    @assert isapprox(basis_inv, reduced_output["inv_state_red"], atol = kry_atol)
    println("reducing model (Krylov) in Julia succeeded")
catch
    println("reducing model (Krylov) in Julia failed")
    Γ1_kry = g1; Γ0_kry = g0
    Π_kry = pi; Ψ_kry = psi
end

println("solving the model...")
T_kry, C_kry, R_kry, ~, ~, ~, ~, eu_kry = gensysct(Γ0_kry, Γ1_kry, C_kry, Ψ_kry, Π_kry,
                                                   complex_decomposition = true)
run_timing_tests && @btime gensysct($Γ0_kry, $Γ1_kry, $C_kry, $Ψ_kry, $Π_kry, complex_decomposition = true);


inflation_kry, monetary_shock_kry, wage_kry, consumption_kry, labor_supply_kry, output_kry =
    simulate_irfs(T_kry, R_kry, vars_SS, periods, steps, Matrix(basis_inv); DisplayLev = DisplayLev)


nothing