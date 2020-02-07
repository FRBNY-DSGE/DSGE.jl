### Model
m = Model904()
TTT, RRR, CCC = solve(m)
meas = measurement(m, TTT, RRR, CCC)
sys = compute_system(m)

nothing
