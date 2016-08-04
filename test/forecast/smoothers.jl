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

pred, vpred, A0, P0 = h5open("$path/../reference/kalman_filter_2part_out.h5", "r") do h5
    read(h5, "pred"), read(h5, "vpred"), read(h5, "z0"), read(h5, "vz0")
end

# Kalman filter test
alpha_hat, eta_hat = kalman_smoother(m, data, TTT, RRR, CCC, QQ, ZZ, DD, A0, P0,
    pred, vpred)

exp_alpha_hat, exp_eta_hat =
    h5open("$path/../reference/kalman_smoother_out.h5", "r") do h5
    read(h5, "alpha_hat"), read(h5, "eta_hat")
end

@test_approx_eq exp_alpha_hat alpha_hat
@test_approx_eq exp_eta_hat eta_hat

# Durbin Koopman smoother test
alpha_hat, eta_hat = durbin_koopman_smoother(m, data, TTT, RRR, CCC, QQ, ZZ, DD,
    A0, P0)

exp_alpha_hat, exp_eta_hat =
    h5open("$path/../reference/durbin_koopman_smoother_out.h5", "r") do h5
    read(h5, "alpha_hat"), read(h5, "eta_hat")
end

@test_approx_eq exp_alpha_hat alpha_hat
@test_approx_eq exp_eta_hat eta_hat


nothing
