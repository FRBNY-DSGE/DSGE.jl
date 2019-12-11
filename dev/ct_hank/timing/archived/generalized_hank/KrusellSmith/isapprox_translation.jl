using JLD, DataStructures, MAT
using ForwardDiff, Plots
using Compat, BenchmarkTools
using BasisMatrices
using DSGE
gr()

### The test file for comparing Julia and Matlab HANK codes in computation of Krusell Smith.


include("../ReductionData_struct.jl")
include("set_parameters.jl")
include("compute_steady_state.jl")
include("equilibrium_conditions.jl")
include("../util.jl")
include("../solve_hjb.jl")
include("../helpers.jl")
include("../simulate.jl")
include("../gensysct.jl")
include("../reduction.jl")

## Step 0: Set Parameters
params, init_params, grids, approx_params, red_params = set_parameters()
atol =  1e-4
println("Parameters set.\n")

## Step 1: Solve for Steady State
println("Computing steady state...")
varsSS = compute_steady_state(grids, params, init_params, approx_params)
vars_SS_mat = matread("saved_outputs/varsSS.mat")

@assert isapprox(vars_SS_mat["VSS"], varsSS[:VSS], atol =  atol)
@assert isapprox(vars_SS_mat["ggSS"], varsSS[:ggSS], atol =  atol)
@assert isapprox(vars_SS_mat["KSS"], varsSS[:KSS], atol =  atol)
@assert isapprox(vars_SS_mat["rSS"], varsSS[:rSS], atol =  atol)
@assert isapprox(vars_SS_mat["wSS"], varsSS[:wSS], atol =  atol)
@assert isapprox(vars_SS_mat["YSS"], varsSS[:YSS], atol =  atol)
@assert isapprox(vars_SS_mat["CSS"], varsSS[:CSS], atol =  atol)
@assert isapprox(vars_SS_mat["ISS"], varsSS[:ISS], atol =  atol)
println("...succeeded!\n")


println("Taking derivatives of equilibrium conditions...")
canonical_form = equilibrium_conditions(varsSS, grids, params, approx_params)
Γ0 = canonical_form[:Γ0]; Γ1 = canonical_form[:Γ1]; C = canonical_form[:C]
Ψ = canonical_form[:Ψ]; Π = canonical_form[:Π]
canonical_preredmat = matread("saved_outputs/canonical_prered.mat")
g0 = canonical_preredmat["g0"]; g1 = canonical_preredmat["g1"]; c = canonical_preredmat["c"]
psi = canonical_preredmat["psi"]; pi = canonical_preredmat["pi"]

try
    @assert isapprox(g0, Γ0, atol =  atol)
    @assert isapprox(c, C, atol =  atol)
    @assert isapprox(psi, Ψ, atol =  atol)
    @assert isapprox(pi, Π, atol =  atol)
    if !isapprox(g1, Γ1, atol =  atol)
        println("Γ1 matrix issue so taking derivatives...")
    end
    @assert isapprox(g1, Γ1, atol = atol)
    println("...succeeded!\n")
catch
    println("...failed.")
    println("Setting canonical matrices to Matlab output.\n")
    Γ1 = full(g1); Γ0 = full(g0); C = full(c); Ψ = full(psi); Π = full(pi)
end

println("Proceeding with reduction steps.")
reduceDistribution = red_params.reduce_distribution
reduceV = red_params.reduce_v

if reduceDistribution == 1
    # Distribution reduction
    println("Reducing distribution using Krylov methods...")
    Γ0_kry, Γ1_kry, Ψ_kry, Π_kry, C_kry, state_red, inv_state_red = krylov_reduction(Γ0, Γ1, Ψ, Π, C, grids, red_params)
    canonical_krymat = matread("saved_outputs/canonical_kry.mat")
    g0 = canonical_krymat["g0"]; g1 = canonical_krymat["g1"]; c = canonical_krymat["c"]
    psi = canonical_krymat["psi"]; pi = canonical_krymat["pi"]
    try
        @assert isapprox(g0, Γ0_kry, atol =  atol)
        @assert isapprox(c, C_kry, atol =  atol)
        @assert isapprox(psi, Ψ_kry, atol =  atol)
        @assert isapprox(pi, Π_kry, atol =  atol)
        if !isapprox(g1, Γ1_kry, atol =  atol)
            println("Γ1 matrix issue so reduction...")
        end
        @assert isapprox(g1, Γ1_kry, atol =  atol)
        println("...succeeded!\n")
        Γ0 = Γ0_kry; Γ1 = Γ1_kry; C = C_kry; Ψ = Ψ_kry; Π = Π_kry
    catch
        println("...failed.")
        println("Setting canonical matrices to Matlab output.\n")
        Γ1 = full(g1); Γ0 = full(g0); C = full(c); Ψ = full(psi); Π = full(pi);
        grids[:n_state_vars_red] = Int64(canonical_krymat["n_g"]); state_red = full(canonical_krymat["state_red"])
        inv_state_red = full(canonical_krymat["inv_state_red"])
    end
else
    grids[:n_state_vars_red] = grids[:n_state_vars]
end

# Need to write a unit test for spline basis
if reduceV == 1
    println("Reducing value function using spline basis...")
    Γ0_spl, Γ1_spl, Ψ_spl, Π_spl, C_spl, to_spline, from_spline = valuef_reduction(grids, red_params, Γ0, Γ1, Ψ, Π, C)
    println("...succeeded!\n")
    Γ0 = Γ0_spl; Γ1 = Γ1_spl; C = C_spl; Ψ = Ψ_spl; Π = Π_spl
elseif reduceV == 0
    from_spline = speye(n_state_vars_red + nVars)
    to_spline = speye(n_state_vars_red + nVars)
    n_splined = nVars
end


## Step 4: Solve Linear System
println("Solving reduced linear system...")
if reduceV && reduceDistribution
    complex_decomp = false
    inverse_basis = full(inv_state_red*from_spline)
elseif reduceDistribution
    complex_decomp = true
    inverse_basis = full(inv_state_red)
else
    complex_decomp = true
end

T, C, R, ~, ~, ~, ~, eu = gensysct(Γ0, Γ1, C, Ψ, Π, complex_decomposition = complex_decomp)

## Step 5: Simulate Impulse Response Functions
println("Simulating IRFs...")

periods = 200
steps = 2000
dt = periods/steps
agg_shock = zeros(1, steps)
agg_shock[1] = 1
sim_states = simulate(T, R, periods, steps, agg_shock; method = :implicit, transformation = inverse_basis)
sim_states = sim_states[grids[:n_jump_vars] + grids[:n_state_vars]:end, :] # extract correct simulated states
sim_states_mat = matread("saved_outputs/sim_states.mat")["simulated"]
inverse_basis_mat = full(matread("saved_outputs/sim_states.mat")["trans_mat"])
post_gensys_mat = matread("saved_outputs/post_gensys.mat")
G1 = full(post_gensys_mat["G1"]); impact = full(post_gensys_mat["impact"])

if isapprox(T, G1, atol = atol) && isapprox(R, impact, atol = atol)
    println("Gensys in Julia succeeded!\n")
else
    println("Gensys in Julia failed.\n")
end

try
    @assert isapprox(sim_states, sim_states_mat, atol =  atol)
    println("IRFs in Julia succeeded using Julia's Gensys outputs!\n")
catch
    println("IRFs in Julia failed using Julia's Gensys outputs.\n")
    try
        sim_states = simulate(G1, impact, periods, steps, agg_shock;
                              method = :implicit, transformation = inverse_basis)
        sim_states = sim_states[grids[:n_jump_vars] + grids[:n_state_vars]:end, :] # extract correct simulated states
        @assert isapprox(sim_states, sim_states_mat, atol =  atol)
        println("IRFs in Julia succeeded using Matlab's Gensys outputs!\n")
    catch
        println("IRFs in Julia failed using Matlab's Gensys outputs.\n")
        try
            sim_states = simulate(G1, impact, periods, steps, agg_shock;
                              method = :implicit, transformation = inverse_basis_mat)
            sim_states = sim_states[grids[:n_jump_vars] + grids[:n_state_vars]:end, :] # extract correct simulated states
            @assert isapprox(sim_states, sim_states_mat, atol =  atol)
            println("IRFs in Julia succeeded using Matlab's Gensys outputs and basis transformation!\n")
        catch
            println("IRFs in Julia failed using Matlab's Gensys outputs and basis transformation.\n")
            sim_states = sim_states_mat
        end
    end
end

println("Calculating impulse response functions...")
irfs_mat = matread("saved_outputs/irfs.mat")
I = length(grids[:a])
tfp = sim_states[1, :] + varsSS[:logAggregateTFP]
output = sim_states[5,:] + varsSS[:YSS]
consumption = sim_states[6,:] + varsSS[:CSS]
investment = sim_states[7,:] + varsSS[:ISS]
tfp_mat = irfs_mat["vAggregateTFP"]'
output_mat = irfs_mat["vAggregateOutput"]'
consumption_mat = irfs_mat["vAggregateConsumption"]'
investment_mat = irfs_mat["vAggregateInvestment"]'

try
    @assert isapprox(tfp, tfp_mat, atol =  atol)
    @assert isapprox(output, output_mat, atol =  atol)
    @assert isapprox(consumption, consumption_mat, atol =  atol)
    @assert isapprox(investment, investment_mat, atol =  atol)
    println("...succeeded!\n")
catch
    println("...failed.\n")
end