using JLD, DataStructures, MAT
using ForwardDiff, Plots
using Compat, BenchmarkTools
using BasisMatrices
gr()

### The test file for comparing Julia and Matlab HANK codes

include("set_parameters.jl")
include("../util.jl")
include("../helpers.jl")
include("compute_steady_state.jl")
include("equilibrium_conditions.jl")
include("../solve.jl")
include("../simulate.jl")
include("../simulate_irfs.jl")
include("../gensysct.jl")
include("../reduction.jl")
include("../valuef_reduction.jl")
include("../ReductionData_struct.jl")
include("../generalized_solver.jl")

# Define dimension size of state space and number of shocks, relations, etc.
I = 100; J = 2
n_v = I*J + 1
n_g = I*J
n_p = 5
n_shocks = 1

params, init_params, grids, approx_params, red_params = set_parameters(n_v, n_g, n_p, n_shocks)
println("Parameters set.")

println("computing steady state...")
vars_SS = compute_steady_state(grids, params, init_params, approx_params)

println("computing canonical form")
canonical_form = equilibrium_conditions(vars_SS, grids, params, approx_params)

#####################
println("computing rest of solution")
generalized_solver(canonical_form, params, init_params, grids, approx_params, red_params)