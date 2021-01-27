using DSGE, JLD2, DifferentialEquations

# Get transition equations
m = KrusellSmithCT()
T, R, C, inverse_basis, basis = solvect(m)

# Make initial guess
n_vars = n_states(m)
n_expectation_errors = n_shocks_expectational(m)
n_shocks = n_shocks_exogenous(m)
x = zeros(Float64, n_vars)

# Set up SDE
f(u,p,t) = T*u
g(u,p,t) = R
dt = 1. / 90 # daily frequency when one period is one month
tspan = (0.0,4. * 50.) # 50 years of data
u0 = reshape(basis * x, size(inverse_basis, 2), 1)
W = WienerProcess(0.0, 0., 0.)
prob = SDEProblem(f, g, u0, tspan, W)
sol = solve(prob, EM(), dt = dt)

save("fine_data_red_basis.jld2", "states", sol.u)
states = sol.u
data = zeros(1, (length(states) - 1) / 90)

# created 2 years worth of daily data -> transform into quarterly data for discrete observations
for i = 1:length(data)
    full_state = inverse_basis * states[1+ (i-1)*90]
    data[i] = full_state[m.endogenous_states[:output][1]]
end
save("simulated_data.jld2", "data", data)
