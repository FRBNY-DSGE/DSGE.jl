using DSGE, JLD

# Solve model
m = KrusellSmithCT()
T, R, C, inverse_basis, basis = solve(m)
T = full(T)
R = full(R)
C = full(C)
inverse_basis = full(inverse_basis)

# Get quarterly data using the truth to change from reduced to full basis
trials = load("many_fine_data_red_basis.jld", "trials")
ntrials = length(keys(trials))
states = trials[1]
data = zeros(1, (length(states) - 1) / 90)

# created daily data -> transform into quarterly data for discrete observations
for i = 1:length(data)
    full_state = inverse_basis * states[1+ (i-1)*90]
    data[i] = full_state[m.endogenous_states[:output][1]]
end

# Initialize storage matrices
σ_vec = collect(linspace(0.005, .009, 20))
σ_vec = sort([σ_vec; .007])
sim_freqs = [2; 3; 12]
loglik_mat = zeros(length(σ_vec), length(sim_freqs), ntrials)
max_σ_mat = zeros(ntrials, length(sim_freqs))
max_lik_mat = zeros(ntrials, length(sim_freqs))
true_lik_vec = zeros(length(sim_freqs))
abs_dist_mat = zeros(ntrials, length(sim_freqs))
st_dev_vec = zeros(length(sim_freqs))
max_abs_dist_vec = zeros(length(sim_freqs))

# Compute true logliik and compile code
for i = 1:length(sim_freqs)
    TT = eye(T) + T * 1/get_setting(m, :state_simulation_freq) # adjust transition equation to account for discretization
    TT, RR, CC = transform_transition_matrices(m, full(TT), full(R), full(C); track_lag = false) # get expanded filter

    # Define measurement equation
    measure = measurement(m, T, R, C, inverse_basis)
    ZZ = measure.ZZ
    DD = measure.DD
    QQ = measure.QQ
    EE = measure.EE

    out = DSGE.kalman_filter(data, TT, RR, CC, QQ, ZZ, DD, EE)
    true_lik_vec[i] = sum(out[1])
end

for j = 1:length(sim_freqs)
    println("Sim_freq: ", sim_freqs[j], " states per quarter")
    update!(m.settings[:state_simulation_freq], Setting(:state_simulation_freq, sim_freqs[j]))

    for n = 1:length(keys(trials))
        states = trials[n]
        data = zeros(1, (length(states) - 1) / 90)

        for i = 1:length(data)
            full_state = inverse_basis * states[1+ (i-1)*90]
            data[i] = full_state[m.endogenous_states[:output][1]]
        end

        for i = 1:length(σ_vec)
            # Recompute steady state solution
            m[:σ_tfp] = σ_vec[i]
            steadystate!(m)
            T, R, C, inverse_basis, basis = solve(m)
            inverse_basis = full(inverse_basis)
            TT = eye(T) + T * 1/get_setting(m, :state_simulation_freq) # adjust transition equation to account for discretization
            TT, RR, CC = transform_transition_matrices(m, TT, R, C; track_lag = false) # get expanded filter

            # Define measurement equation
            measure = measurement(m, T, R, C, inverse_basis)
            ZZ = measure.ZZ
            DD = measure.DD
            QQ = measure.QQ
            EE = measure.EE

            out = DSGE.kalman_filter(data, TT, RR, CC, QQ, ZZ, DD, EE)
            loglik_mat[i, j, n] = sum(out[1])
        end
        max_ind = find(maximum(loglik_mat[:, j, n]) .== loglik_mat[:, j, n])[1]
        max_σ_mat[n, j] = σ_vec[max_ind]
        max_lik_mat[n, j] = maximum(loglik_mat[:, j, n])
    end
    abs_dist_mat[:, j] = max_σ_mat[:, j] - .007
    st_dev_vec[j] = std(abs_dist_mat[:, j])
    max_abs_dist_vec[j] = maximum(abs.(abs_dist_mat[:, j]))

    println("Mean absolute distance")
    println(mean(abs_dist_mat[:, j]))
    println("\n")

    println("Standard deviation of error")
    println(st_dev_vec[j])
    println("\n")

    println("Max absolute error")
    println(max_abs_dist_vec[j])
    println("\n")
end

save("sim_results.jld", "abs_dist_mat", abs_dist_mat, "st_dev_vec", st_dev_vec, "max_abs_dist_vec", max_abs_dist_vec, "max_σ_mat", max_σ_mat, "true_σ", .007, "σ_vec", σ_vec, "max_lik_mat", max_lik_mat, "true_lik_vec", true_lik_vec)
