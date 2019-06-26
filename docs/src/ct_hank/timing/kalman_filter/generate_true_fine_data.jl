using EHANK
using JLD
using DifferentialEquations

# Get transition equations
m = KrusellSmith()
T, R, C, inverse_basis, basis = solve(m)

# Make initial guess
n_vars = get_setting(m, :n_vars)
n_expectation_errors = get_setting(m, :n_expectation_errors)
n_shocks = get_setting(m, :n_shocks)
x = zeros(Float64, n_vars)

# Set up SDE
f(u,p,t) = T*u
g(u,p,t) = R
dt = 1./90 # daily frequency when one period is one quarter
tspan = (0.0,200.) # 50 years of data
u0 = reshape(basis * x, size(inverse_basis, 2), 1)
W = WienerProcess(0.0, 0., 0.)
prob = SDEProblem(f, g, u0, tspan, W)

# Simulate many trials
trials = Dict()
ntrials = 100
for n = 1:ntrials
    sol = solve(prob, EM(), dt = dt)
    trials[n] = sol.u
end

if ntrials > 1
    save("many_fine_data_red_basis.jld", "trials", trials)
else
    save("fine_data_red_basis.jld", "states", trials[1])
end


