using JLD2
m = OneAssetHANK()#=
m <= Setting(:reduce_state_vars, false)
m <= Setting(:reduce_v, false)
=#
T, R, ~, inverse_basis = solvect(m)

n_vars = get_setting(m, :n_vars)
@show n_vars

#TODO: n_expectational_errors does not exist.
#n_expectation_errors = get_setting(m, :n_expectation_errors)
#n_shocks = get_setting(m, :n_shocks)
x = zeros(Float64, n_vars)

nst = size(T,1) ###
#C = zeros(nst)

n = 40 #100
s = zeros(nst,n) #s = zeros(n,406)
s[:,1] = inverse_basis' * x
e = randn(n)

dt = 1/2

for i = 2:n
    s[:,i] = s[:,i-1] + T*s[:,i-1]*dt + R*e[i]*(sqrt(dt))
end
data = zeros(1,Int(n/2))
data[1,:] = [s[64,2*i] for i = 1:Int(n/2)]
save("simulated_data.jld2", "data", data)
