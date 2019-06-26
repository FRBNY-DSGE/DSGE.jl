using EHANK, JLD, BenchmarkTools

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
sim_freqs = [2; 3; 12]
for j = 1:length(sim_freqs)
    println("Sim_freq: ", sim_freqs[j], " states per quarter")
    update!(m.settings[:state_simulation_freq], Setting(:state_simulation_freq, sim_freqs[j]))

    # Run once to compile
    if j == 1
        TT = eye(T) + T * 1/get_setting(m, :state_simulation_freq) # adjust transition equation to account for discretization
        TT, RR, CC = transform_transition_matrices(m, TT, R, C; track_lag = false) # get expanded filter

        # Define measurement equation
        measure = measurement(m, TT, RR, CC, inverse_basis)
        ZZ = measure.ZZ
        DD = measure.DD
        QQ = measure.QQ
        EE = measure.EE

        out = HANK.kalman_filter(data, TT, RR, CC, QQ, ZZ, DD, EE)
    end

    # Define measurement equation
    TT = eye(T) + T * 1/get_setting(m, :state_simulation_freq) # adjust transition equation to account for discretization
    TT, RR, CC = transform_transition_matrices(m, TT, R, C; track_lag = false) # get expanded matrices

    measure = measurement(m, TT, RR, CC, inverse_basis)
    ZZ = measure.ZZ
    DD = measure.DD
    QQ = measure.QQ
    EE = measure.EE

    @btime begin
        ~ = HANK.kalman_filter($data, $TT, $RR, $CC, $QQ, $ZZ, $DD, $EE);
        print()
    end
end
