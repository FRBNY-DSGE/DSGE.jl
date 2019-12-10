using EHANK, JLD, DifferentialEquations, BenchmarkTools
include("../../src/estimate/kalman_filter.jl")

# Solve model
m = KrusellSmith()
T, R, C, inverse_basis, basis = solve(m)

# Get quarterly data using the truth to change from reduced to full basis
states = load("fine_data_red_basis.jld", "states")
data = zeros(1, (length(states) - 1) / 90)

# created daily data -> transform into quarterly data for discrete observations
for i = 1:length(data)
    full_state = inverse_basis * states[1+ (i-1)*90]
    data[i] = full_state[m.endogenous_states[:output][1]]
end

# Compute logliks
σ_vec = collect(linspace(0.001, .013, 20))
σ_vec = sort([σ_vec; .007])
loglik_vec = zeros(length(σ_vec))

# Run once to compile
# Define measurement equation
Δ_t = 1. # length of time traversed b/n observations
dt = 1/.90 # discretization fineness in simulated data
D = [0.]
Z = zeros(1, size(inverse_basis, 1))
Z[m.endogenous_states[:output][1]] = 1
Z = Z * inverse_basis # change to reduced basis
E = zeros(1, 1) * Δ_t # Shock is Browian motion -> variance is Δ_t b/c this is quarterly observation data
Q = eye(1, 1) * dt    # Shock is Brownian motion -> variance is dt b/c this is state data

out = ct_kalman_filter(data, T, R, C, Q, Z, D, E, 1.; method = Euler())


# Now we time Euler and Tsitouras methods
println("Euler method for ODE integration")
@btime begin
~ = ct_kalman_filter($data, $T, $R, $C, $Q, $Z, $D, $E, 1.; method = Euler());
print()
end

println("Tsitouras Runge-Kutta 5/4 method for ODE integration")
@btime begin
~ = ct_kalman_filter($data, $T, $R, $C, $Q, $Z, $D, $E, 1.; method = Tsit5());
print()
end

