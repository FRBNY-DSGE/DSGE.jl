using JLD, DataStructures, MAT
using ForwardDiff, Plots
using Compat, BenchmarkTools
using BasisMatrices
gr()

## Solves the Krusell and Smith (1998)

## Setup the toolbox
# Just need to include folders containing the files in the path.
include("../ReductionData_struct.jl")
include("set_parameters.jl")
include("compute_steady_state.jl")
include("equilibrium_conditions.jl")
include("../generalized_solver.jl")
include("../util.jl")
include("../helpers.jl")
include("../solve.jl")
include("../simulate.jl")
include("../simulate_irfs.jl")
include("../gensysct.jl")
include("../reduction.jl")
include("../generalized_solver.jl")

## Step 0: Set Parameters
params, init_params, grids, approx_params, red_params = set_parameters()
println("Parameters set.")

## Step 1: Solve for Steady State
println("Computing steady state...")
varsSS = compute_steady_state(grids, params, init_params, approx_params)
println("Steady state computed.")

## Step 2: Linearize Model Equations
println("Taking derivatives of equilibrium conditions...")
canonical_form = equilibrium_conditions(varsSS, grids, params, approx_params)
println("Equilibrium conditions computed.")

#println("computing rest of solution")
#generalized_solver(canonical_form, params, approx_params, init_params, grids)

## Step 3. Reducing model
println("Model Reduction ...")

# *State space reduction using Krylov subspace method*
reduceDistribution = red_params.reduce_distribution
reduceV = red_params.reduce_v
Γ1 = canonical_form[:Γ1]; Γ0 = canonical_form[:Γ0]; Π = canonical_form[:Π];
Ψ = canonical_form[:Ψ]; C = canonical_form[:C]

if reduceDistribution == 1
    # State space reduction
    Γ0, Γ1, Ψ, Π, C, state_red,inv_state_red = krylov_reduction(Γ0,Γ1,Ψ, Π, C, grids, red_params)
end

# *Value function reduction using spline inspired bases*
nVars = grids[:nVars]
if reduceV == 1
    # Function calls to create basis reduction
    Γ0, Γ1, Ψ, Π, C, to_spline, from_spline = valuef_reduction(grids, red_params, Γ0, Γ1, Ψ, Π, C)
elseif reduceV == 0
    from_spline = speye(n_g_red + nVars)
    to_spline = speye(n_g_red + nVars)
    n_splined = nVars
end
println("Model reduction completed.")

## Step 4: Solve Linear System
println("Solving reduced linear system...")
if reduceV && reduceDistribution
    complex_decomp = false
    inverse_basis = inv_state_red*from_spline
elseif reduceDistribution
    complex_decomp = true
    inverse_basis = inv_state_red
else
    complex_decomp = true
end

T, C, R, ~, ~, ~, ~, eu = gensysct(Γ0, Γ1, C, Ψ, Π, complex_decomposition = complex_decomp)
println("Solving reduced linear system completed.")

## Step 5: Simulate Impulse Response Functions
println("Simulating Model...")

# initialize shocks for simulation
periods = 200
steps = 2000
shocks = Matrix(inverse_basis)
global shocks, R, T, periods, steps

simulate(T, R, periods, steps, Matrix(inverse_basis); method=:implicit)
println("IRFs calculated.")


