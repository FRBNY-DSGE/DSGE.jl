using JLD2
using DSGE
m = KrusellSmithCT()
T, R, ~, inverse_basis = solvect(m)

n_vars = get_setting(m, :n_vars)
n_expectation_errors = get_setting(m, :n_expectation_errors)
n_shocks = get_setting(m, :n_shocks)
x = zeros(Float64, n_vars)

nst = size(T,1) ###
C = zeros(nst)

n = 40 #100
s = zeros(nst,n) #s = zeros(n,406)
s[:,1] = inverse_basis' * x
e = randn(n)

dt = 1/2

for i = 2:n
    s[:, i] = s[:, i-1] + T*s[:,i-1]*dt + R*e[i]*(sqrt(dt))
end
data = zeros(1,Int(n/2))
data[1, :] = [s[nst - 2, 2 * i] for i = 1:Int(n/2)]

save("simulated_data.jld2", "data", data)
