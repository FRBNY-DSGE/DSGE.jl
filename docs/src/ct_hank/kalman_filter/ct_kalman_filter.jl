using EHANK, JLD, DifferentialEquations
include("../../src/estimate/kalman_filter.jl")

# This file will test how well this Kalman Filter finds the
# true outcome when we generate new data many different times

# Solve model
m = KrusellSmith()
T, R, C, inverse_basis, basis = solve(m)

# Get quarterly data using the truth to change from reduced to full basis
trials = load("many_fine_data_red_basis.jld", "trials")
ntrials = length(keys(trials))
states = trials[1]
data = zeros(1, (length(states) - 1) / 90)

# created 50 years worth of daily data -> transform into quarterly data for discrete observations
for i = 1:length(data)
    full_state = inverse_basis * states[1+ (i-1)*90]
    data[i] = full_state[m.endogenous_states[:output][1]]
end

# Compute logliks
σ_vec = collect(linspace(0.005, .009, 20))
σ_vec = sort([σ_vec; .007])
loglik_mat = zeros(length(σ_vec), ntrials)

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

out = ct_kalman_filter(data, T, R, C, Q, Z, D, E, 1.)
true_lik = sum(out[1])

# Start timing
max_σ = zeros(ntrials)
max_lik = zeros(ntrials)
for n = 1:ntrials
    states = trials[n] # get trial data
    data = zeros(1, (length(states) - 1) / 90)
    # created 2 years worth of daily data -> transform into quarterly data for discrete observations
    for i = 1:length(data)
        full_state = inverse_basis * states[1+ (i-1)*90]
        data[i] = full_state[m.endogenous_states[:output][1]]
    end

    for i = 1:length(σ_vec)
        # Recompute steady state solution
        m[:σ_tfp] = σ_vec[i]
        steadystate!(m)
        T, R, C, inverse_basis, basis = solve(m)

        # Define measurement equation
        Δ_t = 1. # length of time traversed b/n observations
        dt = 1/.90 # discretization fineness in simulated data
        D = [0.]
        Z = zeros(1, size(inverse_basis, 1))
        Z[m.endogenous_states[:output][1]] = 1
        Z = Z * inverse_basis # change to reduced basis
        E = zeros(1, 1) * Δ_t # Shock is Browian motion -> variance is Δ_t b/c this is quarterly observation data
        Q = eye(1, 1) * dt    # Shock is Brownian motion -> variance is dt b/c this is state data
        out = ct_kalman_filter(data, T, R, C, Q, Z, D, E, 1.)
        loglik_mat[i, n] = sum(out[1])
    end
    max_ind = find(maximum(loglik_mat[:, n]) .== loglik_mat[:, n])[1]
    max_lik[n] = maximum(loglik_mat[:, n])
    max_σ[n] = σ_vec[max_ind]
end

# Report mean abs distance (true lik sigma- estimated max lik sigma), st dev, and max
abs_dist = max_σ - .007
st_dev = std(max_σ - .007)
max_abs_dist = maximum(abs.(max_σ - .007))

println("Mean absolute distance")
println(mean(abs_dist))

println("\n")

println("Standard deviation of error")
println(st_dev)
println("\n")

println("Max absolute error")
println(max_abs_dist)

save("ct_results.jld", "abs_dist", abs_dist, "st_dev", st_dev, "max_abs_dist", max_abs_dist, "max_σ", max_σ, "true_σ", .007, "σ_vec", σ_vec, "max_lik", max_lik, "true_lik", true_lik)
