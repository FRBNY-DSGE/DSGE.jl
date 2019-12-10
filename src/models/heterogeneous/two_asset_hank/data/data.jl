#Simulate data series for inflation
using HDF5
using StateSpaceRoutines

m = OneAssetHANK()

n_vars = get_setting(m, :n_vars)
n_expectation_errors = get_setting(m, :n_expectation_errors)
n_shocks = get_setting(m, :n_shocks)
x = zeros(Float64, n_vars)

I = get_setting(m, :I)
n_v = get_setting(m, :n_jump_vars)
n_g = get_setting(m, :n_state_vars)

x[1:n_v-1] = reshape(m[:V_ss].value, 2*I, 1)
x[n_v] = m[:inflation_ss].value
x[n_v + 1 : n_v + n_g - 1] = m[:g_ss].value[1:end-1]
x[n_v + n_g + 1] = m[:w_ss].value
x[n_v + n_g + 2] = m[:N_ss].value
x[n_v + n_g + 3] = m[:C_ss].value
x[n_v + n_g + 4] = m[:Y_ss].value
x[n_v + n_g + 5] = m[:B_ss].value

T, R, ~, inverse_basis = solvect(m)

nst = size(T,1) ###
C = zeros(nst)
TTT, RRR, C = transform_transition_matrices(m, T, R, C, true)

n = 40#100
s = zeros(nst,n) #s = zeros(n,406)
s[:,1] = inverse_basis'*x
e = randn(n)

dt = 1/2

for i = 2:n
  s[:,i] = T*s[:,i-1]*dt + R*e[i]*(sqrt(dt))
end

data = zeros(1,Int(n/2))
data[1,:] = [s[64,2*i] for i = 1:Int(n/2)]

#Z is 1x130

#write("output.h5", "data", data)

#Create error covariance matrix, matrices in measurement equation


Q = eye(2,2)#
E = zeros(1,1)#
   #steady state output
D = [mean(data)]
Z = zeros(1, nst*2)#(1, nst) # to allow for reduction * inverse_basis #
   #where inflation series shows up in matrix of states
Z[nst*2-1] = 1

s_0 = zeros(nst*2)
   #Exceeds max iterations if we use discrete lyapunov solver
P_0 = 1e6*eye(nst*2) #solve_discrete_lyapunov(T_tilde, R*Q*R')

#@assert false

#Call standard kalman filter
output = kalman_filter(data, TTT, RRR, C, Q, Z, D, E, s_0, P_0)


###
#=

=#