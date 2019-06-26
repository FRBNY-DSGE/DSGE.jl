using JLD, HDF5
using QuantEcon: solve_discrete_lyapunov
using StateSpaceRoutines

#Grab transition equation matrices from model solution
state_space_vars = load("data/state_space_vars.jld")
T = state_space_vars["T"]
R = state_space_vars["R"]
C = vec(state_space_vars["C"])

#For weekly observations
convert = 13
dt = 1/convert

#grab inflation series
data_orig = h5read("data/smc.h5","data")[2,:]

#Size of original quarterly series
nt = length(data_orig)
#Size of series with weekly observations
Nt = convert*nt

#Create error covariance matrix, matrices in measurement equation
Q = eye(1,1)*dt
E = zeros(1,1)
   #steady state inflation
D = [0.0]
Z = zeros(1, 401) # to allow for reduction * inverse_basis
   #where inflation series shows up in matrix of states
Z[201] = 400

#Create expanded matrix of data to accomodate weekly observations
data = fill(NaN, 1, Nt)
for i in 1:Nt
    if i%13 == 1
        data[i] = data_orig[Int(ceil(i/13))]
    end
end

#Create new transition matrix: (T*Δt + I)
T_tilde = T*dt + eye(size(T,1))

n_states = size(T, 1)
s_0 = zeros(n_states)
   #Exceeds max iterations if we use discrete lyapunov solver
P_0 = 1e6*eye(n_states) #solve_discrete_lyapunov(T_tilde, R*Q*R')

@assert false

#Call standard calman filter
output = kalman_filter(data, T_tilde, R, C, Q, Z, D, E, s_0, P_0)
#Print loglikelihood
println(sum(output[1]))

