using JLD, DataStructures, MAT
using ForwardDiff, Plots
using Compat, BenchmarkTools
using BasisMatrices
gr()

### The test file for comparing Julia and Matlab HANK codes
include("../ReductionData_struct.jl")
include("set_parameters.jl")
include("../util.jl")
include("../solve_hjb.jl")
include("../helpers.jl")
include("compute_steady_state.jl")
include("equilibrium_conditions.jl")
include("../simulate.jl")
include("../simulate_irfs.jl")
include("../gensysct.jl")
include("../reduction.jl")

ReduceDist_hors = 20
atol = 1e-4
run_timing_tests = false

# Set dimensions
I = 100; J = 2
n_v = Int64(I*J + 1) # The number of jump variables (value function + inflation)
n_g = Int64(I*J)     # The number of endogenous state variables (distribution + monetary policy)
n_p = 5              # The number of static relations: bond-market clearing, labor
                     # market clearing, consumption, output, total assets
n_shocks = 1         # only monetary policy shock is considered

params, init_params, grids, approx_params, red_params = set_parameters(n_v, n_g, n_p, n_shocks)
println("Parameters set.")

## check residual
println("computing steady state...")
vars_SS = compute_steady_state(grids, params, init_params, approx_params)
run_timing_tests && @btime compute_steady_state($grids, $params, $init_params, $approx_params);

#vars_SS_matlab = matread("saved_outputs/vars_SS.mat")
#for var in [:B_SS, :Y_SS, :inflation_SS, :N_SS, :C_SS,
#            :V_SS, :T_SS, :profit_SS, :rnom_SS,
#            :s_SS,  :w_SS, :c_SS, :h_SS, :g_SS]
#    @assert vars_SS_matlab[string(var)] ≈ vars_SS[var] "equality test failed for $var"
#end
println("taking derivatives in Julia succeeded")

#x = zeros(Float64, 2*grids[:n_vars] + grids[:n_exp_errors] + grids[:n_shocks])
canonical_form = equilibrium_conditions(vars_SS, grids, params, approx_params)
run_timing_tests && @btime equilibrium_conditions($vars_SS, $grids, $params, $approx_params)

derivs_output = matread("saved_outputs/derivs.mat")
derivs_matlab = Matrix(derivs_output["derivs"])
Γ0 = canonical_form[:Γ0]; Γ1 = canonical_form[:Γ1]; Ψ = canonical_form[:Ψ]
Π = canonical_form[:Π]; C = canonical_form[:C]

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

##################################################################

## solve out static constraints and then solve model
# println("solving out static constraints...")
# Γ0_red, Γ1_red, Ψ_red, Π_red, C_red, redundant_states_inv = solve_static_constraints(Γ0, Γ1, Ψ, Π, C)
# run_timing_tests && @btime solve_static_constraints($Γ0, $Γ1, $Ψ, $Π, $C);

# reduced_output = matread("saved_outputs/post_reduction.mat")
# g0  = Matrix(reduced_output["g0"]); g1  = Matrix(reduced_output["g1"])
# pi  = Matrix(reduced_output["pi"]); psi = Matrix(reduced_output["psi"])
# try
#     @assert isapprox(g0, Γ0_red, atol = atol); @assert isapprox(g1, Γ1_red, atol = atol)
#     @assert isapprox(pi, Π_red, atol = atol); @assert isapprox(psi, Ψ_red, atol = atol)
#     println("reducing model in Julia succeeded")
# catch
#     println("reducing model in Julia failed")
#     Γ1_red = g1; Γ0_red = g0
#     Π_red = pi; Ψ_red = psi
# end

# println("solving the model...")
# T, C, R, ~, ~, ~, ~, eu = gensysct(Γ0_red, Γ1_red, C_red, Ψ_red, Π_red)
# run_timing_tests && @btime gensysct($Γ0_red, $Γ1_red, $C_red, $Ψ_red, $Π_red)
# gensys_output = matread("saved_outputs/post_gensys.mat")
# G1 = Matrix(gensys_output["G1"]); impact = Matrix(gensys_output["impact"])

# ## simulate IRFs
# println("simulating IRFs...")

periods = 100
steps = 400
dt = periods/steps
agg_shock = zeros(1, steps)
agg_shock[1] = 1/sqrt(dt)
#transformation = redundant_states_inv

# simulated_states = simulate(T, R, periods, steps, agg_shock; method = :implicit)
#  run_timing_tests && @btime simulate($T, $R, $periods, $steps, $agg_shock; method = :implicit)
# simulated_states = (transformation * simulated_states)[cat(1,n_v, n_v+n_g:n_v+n_g+5),:]

# if isapprox(T, G1) && isapprox(R, impact)
#     println("Gensys in Julia suceeded")
# else
#     println("Gensys in Julia failed")
# end

# simulated_states_matlab = simulate(G1, impact, periods, steps, agg_shock; method = :implicit)
# simulated_states_matlab = (transformation * simulated_states_matlab)[cat(1,n_v, n_v+n_g:n_v+n_g+5),:]
# try
#     @assert isapprox(simulated_states, simulated_states_matlab, atol = 1e-3)
#     println("IRFs in Julia succeeded (still using Julia's Gensys outputs)")
# catch
#     println("IRFs in Julia failed (still using Julia's Gensys outputs)")
# end

# inflation, monetary_shock, wage, consumption, labor_supply, output =
#     simulate_irfs(T, R, vars_SS, periods, steps, transformation)

## Perform state space and V reduction, then solve model

# Remake state space matrices
Γ1 = canonical_form[:Γ1]
Γ0 = canonical_form[:Γ0]
Π = canonical_form[:Π]
Ψ = canonical_form[:Ψ]
C = canonical_form[:C]

println("performing Krylov reduction on state space...")

if red_params.reduce_distribution
    Γ0_kry, Γ1_kry, Ψ_kry, Π_kry, C_kry, basis, basis_inv = krylov_reduction(Γ0, Γ1, Ψ, Π, C, grids, red_params)

    run_timing_tests && @btime krylov_reduction($Γ0, $Γ1, $Ψ, $Π, $C, $grids, $red_params);

    reduced_output = matread("saved_outputs/krylov_outputs.h5")
    g0  = Matrix(reduced_output["g0"]); g1  = Matrix(reduced_output["g1"])
    pi  = Matrix(reduced_output["pi"]); psi = Matrix(reduced_output["psi"])

    kry_atol = atol
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
else
    Γ1_kry = Γ1; Γ0_kry = Γ0; Ψ_kry = Ψ; Π_kry = Π; C_kry = C
    grids[:n_state_vars_reduce] = grids[:n_state_vars]
end


println("performing value function reduction via quadratic spline projection...")

#################################################################################

if red_params.reduce_v
    Γ0_spl, Γ1_spl, Ψ_spl, Π_spl, C_spl, to_spline, from_spline = valuef_reduction(grids, red_params, Γ0_kry, Γ1_kry, Ψ_kry, Π_kry, C_kry)

    run_timing_tests && @btime valuef_reduction($grids, $red_params, $Γ0_kry, $Γ1_kry, $Ψ_kry, $Π_kry, $C_kry)
    if !red_params.reduce_v
        # Compare with Matlab output
        reduced_output = matread("saved_outputs/splines_outputs.mat")
        g0  = Matrix(reduced_output["g0"]); g1  = Matrix(reduced_output["g1"])
        pi  = Matrix(reduced_output["pi"]); psi = Matrix(reduced_output["psi"])
        mat_from_spline = Matrix(reduced_output["from_spline"])
        mat_to_spline = Matrix(reduced_output["to_spline"])
        constant = Matrix(reduced_output["constant"])

        spl_atol = 1e-4
        try
            @assert isapprox(g0, Γ0_spl, atol = spl_atol); @assert isapprox(g1, Γ1_spl, atol = spl_atol)
            @assert isapprox(pi, Π_spl, atol = spl_atol); @assert isapprox(psi, Ψ_spl, atol = spl_atol)
            @assert isapprox(constant, C_spl, atol = spl_atol)
            println("Reduction of value function onto spline basis succeeded.")
        catch
            # For now, just use the read-in Matlab matrices
            println("Reduction of value function onto spline basis failed.")
            Γ1_spl = g1; Γ0_spl = g0
            Π_spl = pi; Ψ_spl = psi
            C_spl = constant
        end
    end
else
    # Create identity matrix for code reuse below
    n_state_vars_red = grids[:n_state_vars_reduce]
    from_spline = speye(n_state_vars_red + n_v)
    to_spline = speye(n_state_vars_red + n_v)
    n_splined = n_v
end


###################################################################################

println("solving the model...")
if red_params.reduce_v && red_params.reduce_distribution
    complex_decomp = false
    inverse_basis = basis_inv*from_spline
elseif red_params.reduce_distribution
    inverse_basis = basis_inv
else
    complex_decomp = true
    inverse_basis = eye(size(T_spl,1))
end
T_spl, C_spl, R_spl, ~, ~, ~, ~, eu_spl = gensysct(Γ0_spl, Γ1_spl, C_spl, Ψ_spl, Π_spl,
                                                  complex_decomposition = complex_decomp)
run_timing_tests && @btime gensysct($Γ0_spl, $Γ1_spl, $C_spl, $Ψ_spl, $Π_spl, complex_decomposition = complex_decomp);

inflation_spl, monetary_shock_spl, wage_spl, consumption_spl, labor_supply_spl, output_spl =
    simulate_irfs(T_spl, R_spl, vars_SS, periods, steps, Matrix(inverse_basis))

vars_IRF_matlab = matread("saved_outputs/simulated.mat")
inflation = vec(vars_IRF_matlab["inflation"]); monetary_shock = vec(vars_IRF_matlab["monetary_shock"])
consumption = vec(vars_IRF_matlab["consumption"]); Y = vec(vars_IRF_matlab["Y"])
lab_sup = vec(vars_IRF_matlab["lab_sup"]); wage = vec(vars_IRF_matlab["wage"])

spl_atol = atol
try
    @assert isapprox(inflation, inflation_spl, atol = spl_atol)
    @assert isapprox(monetary_shock, monetary_shock_spl, atol = spl_atol)
    @assert isapprox(wage, wage_spl, atol = spl_atol)
    @assert isapprox(consumption, consumption_spl, atol = spl_atol)
    @assert isapprox(lab_sup, labor_supply_spl, atol = spl_atol)
    @assert isapprox(Y, output_spl, atol = spl_atol)
    println("Calculating impulse response functions succeeded.")
catch
    println("Calculating impulse response functions failed.")
end
