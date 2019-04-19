using JLD
using DSGE
using DifferentialEquations
using ForwardDiff

# Get transition equations
m = KrusellSmithCT()
T, R, C, inverse_basis, basis = DSGE.solve(m)

# Make initial guess
n_vars = n_states(m)
n_expectation_errors = n_shocks_expectational(m)
n_shocks = n_shocks_exogenous(m)
x = zeros(Float64, n_vars)

# Set up SDE
f(u,p,t) = T*u
g(u,p,t) = R
dt = 1 ./ 90 # daily frequency when one period is one quarter
tspan = (0.0,200.) # 50 years of data
u0 = reshape(basis * x, (size(inverse_basis, 2), 1))
W = WienerProcess(0.0, 0., 0.)
prob = SDEProblem(f, g, u0, tspan, W)

# Simulate many trials
trials = Dict()
ntrials = 100
for n = 1:ntrials
    sol = DifferentialEquations.solve(prob, EM(), dt = dt)
    trials[n] = sol.u
end

if ntrials > 1
    save("many_fine_data_red_basis.jld", "trials", trials)
else
    save("fine_data_red_basis.jld", "states", trials[1])
end
