using JLD, DataStructures, MAT
using ForwardDiff, Plots
import Gensys: gensysct
gr()

# The spec file for computing the steady state of the one asset HANK model

include("set_parameters.jl")
include("util.jl")
include("helpers.jl")
include("compute_steady_state.jl")
include("equilibrium_conditions.jl")
include("solve.jl")
include("simulate.jl")
include("simulate_irfs.jl")
include("gensysct.jl")

const DisplayLev = 1

# Solve model
T, C, R, vars_SS, redundant_states_inv = solve(DisplayLev = DisplayLev)

# Simulate IRFs
periods = 100
steps = 400
transformation = redundant_states_inv
inflation, monetary_shock, wage, consumption, labor_supply, output =
    simulate_irfs(T, R, vars_SS, periods, steps, transformation; DisplayLev = DisplayLev)