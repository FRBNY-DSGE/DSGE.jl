using BenchmarkTools
using ClusterManagers
using HDF5, JLD, Base.Test
using DataFrames
using DSGE
import DSGE.update!
using QuantEcon: solve_discrete_lyapunov, solve_discrete_riccati
using StateSpaceRoutines

#include("large_scale_kalman.jl")

m = AnSchorfheide()

data = h5read("../data/smc.h5","data")
params = [2.09, 0.98, 2.25, 0.65, 0.34, 3.16, 0.51, 0.81, 0.98, 0.93, 0.19, 0.65, 0.24,
          0.115985, 0.294166, 0.447587]
update!(m, params)

# Solution to a Linear DSGE Model w/ IID Gaussian Errors

system  = compute_system(m)
T     = system[:TTT]
R     = system[:RRR]
C     = system[:CCC]
Q      = system[:QQ]
Z      = system[:ZZ]
D      = system[:DD]
E      = system[:EE]

#test stability of T (if T [in their notation F] according to http://webee.technion.ac.il/people/shimkin/Estimation09/ch6_ss.pdf THM 1 P approaches Pbar and Pbar is the unique non-negative-definite solution of the Algebraic Riccati Equation) For stability |λ_i| < 1 for all i
if all(eigvals(T) .< 1)
    println("T is stable.")
    #Pbar = solve_discrete_riccati(A, B, R, Q, N) see http://quantecon.github.io/QuantEcon.jl/latest/api/QuantEcon.html and compare inputs to above link
else
    println("T is not stable.")
end

# Generation of the initial state draws
n_states = n_states_augmented(m)
s_0 = zeros(n_states)
P_0 = solve_discrete_lyapunov(T, R*Q*R')
@assert false

#Kalman
#norm_P_T, ch_ll, truelik = compute_values(data, T, R, C, Q, Z, D, E, s_0, P_0)

#kaloutput = Dict{Symbol, Any}()

#kaloutput[:norm_P_T] = norm_P_T
#kaloutput[:ch_ll] = ch_ll
#kaloutput[:truelik] = truelik

# Kalman
# compute_values(data, T, R, C, Q, Z, D, E, s_0, P_0)

# Processes for g_t and z_t are both AR(1) and exogenous;
# monetary policy shock is endogenous b/c depends on endogenous rules,
# despite the AR(1) shock structure
M = eye(n_states)
M[1, 1] = 0; M[1, 6] = 1
M[6, 6] = 0; M[6, 1] = 1
M[2, 2] = 0; M[2, 5] = 1
M[5, 5] = 0; M[5, 2] = 1
# M[1,1] = 0; M[1, 5] = 1
# M[5, 5] = 0; M[5, 1] = 1
# M[2, 2] = 0; M[2, 6] = 1
# M[6, 6] = 0; M[6, 2] = 1
Mtild = eye(3) # We preserve orderings of state shocks

# Zero out entries that should be zero
T_ord = M' * T * M
R_ord = M' * R * Mtild
T_ord[1:2, 3:end] = 0; T_ord[1, 2] = 0; T_ord[2, 1] = 0
R_ord[1, 2] = 0; R_ord[2, 1] = 0
T = inv(M') * T_ord * inv(M)
R = inv(M') * R_ord * inv(Mtild)
output_true = kalman_filter(data, T, R, C, Q, Z, D, E, s_0, P_0)
#output = block_kalman_filter(data, T, R, C, Q, Z, D, E, M, Mtild, [2; 0; 0; 6], s_0, P_0)
println(sum(output_true[1]))
#println(sum(output[1]))

