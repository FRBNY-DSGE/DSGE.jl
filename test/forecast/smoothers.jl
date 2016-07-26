using DSGE, HDF5, Base.Test
include("../util.jl")

path = dirname(@__FILE__)


# Kalman smoother test
h5 = h5open("$path/../reference/kalman_smoother_args.h5")
for arg in ["A0", "C", "P0", "Q", "R", "T", "Z", "antlags", "b", "nant",
            "peachcount", "pred", "psize", "vpred", "y"]
    eval(parse("$arg = read(h5, \"$arg\")"))
end
close(h5)

m = Model990()
m.testing = true
m <= Setting(:n_anticipated_shocks, 6, "Number of anticipated policy shocks")
DSGE.init_model_indices!(m)
m <= Setting(:date_forecast_start, quartertodate("2016-Q1"))
m <= Setting(:n_anticipated_lags, antlags)

smoothed = kalman_smoother(m, y, T, R, C, Q, Z, b, A0, P0, pred, vpred)
alpha_hat, eta_hat = smoothed.states, smoothed.shocks

Nt0 = n_presample_periods(m)
exp_alpha_hat, exp_eta_hat =
    h5open("$path/../reference/kalman_smoother_out.h5", "r") do h5
        read(h5, "alpha_hat")[:, (Nt0+1):end], read(h5, "eta_hat")[:, (Nt0+1):end]
end

@test_approx_eq exp_alpha_hat alpha_hat
@test_approx_eq exp_eta_hat eta_hat


# Durbin Koopman smoother test
data, P0 = h5open("$path/../reference/durbin_koopman_smoother_args.h5") do h5
    read(h5, "YY_all")', read(h5, "P0")
end

m = Model990()
m.testing = true
m <= Setting(:n_anticipated_shocks, 6, "Number of anticipated policy shocks")
DSGE.init_model_indices!(m)
m <= Setting(:date_zlbregime_start, quartertodate("2008-Q3"))
m <= Setting(:date_mainsample_end, quartertodate("2015-Q1"))
m <= Setting(:date_forecast_start, quartertodate("2015-Q2"))
m <= Setting(:n_anticipated_lags, DSGE.subtract_quarters(get_setting(m, :date_forecast_start), get_setting(m, :date_zlbregime_start)))

TTT, RRR, CCC = solve(m)
meas = measurement(m, TTT, RRR, CCC)
QQ, ZZ, DD = meas.QQ, meas.ZZ, meas.DD
A0 = zeros(size(P0, 1))

smoothed = durbin_koopman_smoother(m, data, TTT, RRR, CCC, QQ, ZZ, DD, A0, P0)

Nt0 = n_presample_periods(m)
exp_alpha_hat, exp_eta_hat =
    h5open("$path/../reference/durbin_koopman_smoother_out.h5", "r") do h5
    read(h5, "alpha_hat")[:, Nt0+1:end], read(h5, "eta_hat")[:, Nt0+1:end]
end
alpha_hat, eta_hat = smoothed.states, smoothed.shocks

@test_approx_eq exp_alpha_hat alpha_hat
@test_approx_eq exp_eta_hat eta_hat


nothing
