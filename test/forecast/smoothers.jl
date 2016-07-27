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
m <= Setting(:date_mainsample_end, quartertodate("2015-Q4"))
m <= Setting(:date_forecast_start, quartertodate("2016-Q1"))
m <= Setting(:n_anticipated_lags, DSGE.subtract_quarters(get_setting(m, :date_forecast_start),
    get_setting(m, :date_zlbregime_start)) - 1)
@assert n_anticipated_lags(m) == 28

TTT, RRR, CCC = solve(m)
meas = measurement(m, TTT, RRR, CCC)
QQ, ZZ, DD = meas.QQ, meas.ZZ, meas.DD

A0 = zeros(size(TTT, 1))
P0 = QuantEcon.solve_discrete_lyapunov(TTT, RRR*QQ*RRR')
R2, R3, R1 = kalman_filter_2part(m, data', TTT, RRR, CCC, A0, P0, allout = true, augment_states = true)
pred  = hcat(R1[:pred], R2[:pred], R3[:pred])
vpred = cat(3, R1[:vpred], R2[:vpred], R3[:vpred])


# Kalman filter test
smoothed = kalman_smoother(m, data, TTT, RRR, CCC, QQ, ZZ, DD, A0, P0, pred, vpred)
alpha_hat, eta_hat = smoothed.states, smoothed.shocks

exp_alpha_hat, exp_eta_hat =
    h5open("$path/../reference/kalman_smoother_out.h5", "r") do h5
    read(h5, "alpha_hat"), read(h5, "eta_hat")
end

@test_approx_eq exp_alpha_hat alpha_hat
@test_approx_eq exp_eta_hat eta_hat


# Durbin Koopman smoother test
smoothed = durbin_koopman_smoother(m, data, TTT, RRR, CCC, QQ, ZZ, DD, A0, P0)
alpha_hat, eta_hat = smoothed.states, smoothed.shocks

exp_alpha_hat, exp_eta_hat =
    h5open("$path/../reference/durbin_koopman_smoother_out.h5", "r") do h5
    read(h5, "alpha_hat"), read(h5, "eta_hat")
end

@test_approx_eq exp_alpha_hat alpha_hat
@test_approx_eq exp_eta_hat eta_hat


nothing
