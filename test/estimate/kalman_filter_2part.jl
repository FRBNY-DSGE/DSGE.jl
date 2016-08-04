using DSGE, HDF5, Base.Test
include("../util.jl")

path = dirname(@__FILE__)

# Set up
data = h5open("$path/../reference/smoother_args.h5", "r") do h5
    read(h5, "data")
end

m = Model990()
m.testing = true
m <= Setting(:n_anticipated_shocks, 6)
DSGE.init_model_indices!(m)
m <= Setting(:date_forecast_start, quartertodate("2016-Q1"))

TTT, RRR, CCC = solve(m)
meas = measurement(m, TTT, RRR, CCC)
QQ, ZZ, DD = meas.QQ, meas.ZZ, meas.DD

state_inds = inds_states_no_ant(m)
shock_inds = inds_shocks_no_ant(m)
TTT_small = TTT[state_inds, state_inds]
RRR_small = RRR[state_inds, shock_inds]
QQ_small  = QQ[shock_inds, shock_inds]

A0 = zeros(n_states_augmented(m))
P0_small = QuantEcon.solve_discrete_lyapunov(TTT_small, RRR_small*QQ_small*RRR_small')
P0 = zeros(n_states_augmented(m), n_states_augmented(m))
P0[state_inds, state_inds] = P0_small

# Kalman filter with all arguments provided
k1, _, _, _ = kalman_filter_2part(m, data, TTT, RRR, CCC, A0, P0; allout = true,
    include_presample = true)

# Kalman filter without z0 and vz0
k2, _, _, _ = kalman_filter_2part(m, data, TTT, RRR, CCC; allout = true,
    include_presample = true)

# Test equality
h5 = h5open("$path/../reference/kalman_filter_2part_out.h5", "r")
for out in [:L, :zend, :Pend, :pred, :vpred, :yprederror, :ystdprederror, :rmse,
            :rmsd, :filt, :vfilt, :z0, :vz0]
    expect = read(h5, "$out")

    @test_approx_eq expect k1[out]
    @test_approx_eq expect k2[out]
end
close(h5)

nothing