
using BenchmarkTools
using ClusterManagers
using HDF5, JLD, Base.Test
using DataFrames
using DSGE
import DSGE.update!
using QuantEcon: solve_discrete_lyapunov, solve_discrete_riccati

include("test.jl")
include("kalman_filter.jl")
include("ar1_ct_kalmanfunction.jl")

data, T, R, C, Q, Z, D, E, n_states = get_oneassethank_statespace()

if all(eigvals(T) .< 1)
    println("T is stable.")
else
    println("T is not stable.")
end

# Generation of the initial state draws
s_0 = zeros(n_states)
P_0 = solve_discrete_lyapunov(T, R*Q*R')

#Kalman
norm_P_T, ch_ll, truelik = compute_values(data, T, R, C, Q, Z, D, E, s_0, P_0)

kaloutput = Dict{Symbol, Any}()

kaloutput[:norm_P_T] = norm_P_T
kaloutput[:ch_ll] = ch_ll
kaloutput[:truelik] = truelik