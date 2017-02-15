using DSGE, HDF5, JLD, Base.Test
include("../util.jl")

path = dirname(@__FILE__)

# Set up
data, TTT, RRR, CCC = h5open("$path/../reference/kalman_filter_2part_args.h5", "r") do file
    read(file, "data"), read(file, "TTT"), read(file, "RRR"), read(file, "CCC")
end
kal = jldopen("$path/../reference/kalman_filter_2part_out.jld", "r") do file
    read(file, "exp_kal")
end

custom_settings = Dict{Symbol, Setting}(
    :n_anticipated_shocks => Setting(:n_anticipated_shocks, 6),
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2016-Q1")))
m = Model990(custom_settings = custom_settings, testing = true)

meas = measurement(m, TTT, RRR, CCC)
QQ, ZZ, DD, MM, EE = meas.QQ, meas.ZZ, meas.DD, meas.MM, meas.EE


# Read expected output
exp_alpha_hat, exp_eta_hat = h5open("$path/../reference/kalman_smoother_out.h5", "r") do h5
    read(h5, "alpha_hat"), read(h5, "eta_hat")
end


# Kalman smoother with anticipated shocks
alpha_hat, eta_hat = kalman_smoother(m, data, TTT, RRR, CCC, QQ, ZZ, DD, kal[:z0], kal[:vz0],
    kal[:pred], kal[:vpred])

@test_approx_eq exp_alpha_hat alpha_hat
@test_approx_eq exp_eta_hat eta_hat
@test_approx_eq kal[:zend] alpha_hat[:, end]

# Hamilton smoother with anticipated shocks
system = compute_system(m)
alpha_hat, eta_hat = hamilton_smoother(m, data, system, kal)

exp_alpha_hat, exp_eta_hat = h5open("$path/../reference/Hamilton_smoother_out.h5", "r") do h5
    read(h5, "alpha_hat"), read(h5, "eta_hat")
end

@test_approx_eq exp_alpha_hat alpha_hat
@test_approx_eq exp_eta_hat eta_hat
@test_approx_eq kal[:zend] alpha_hat[:, end]


# Carter and Kohn smoother with anticipated shocks
system = compute_system(m)
alpha_hat, eta_hat = carter_kohn_smoother(m, data, system, kal)


exp_alpha_hat, exp_eta_hat = h5open("$path/../reference/Hamilton_smoother_out.h5", "r") do h5
    read(h5, "alpha_hat"), read(h5, "eta_hat")
end

@test_approx_eq exp_alpha_hat alpha_hat
@test_approx_eq exp_eta_hat eta_hat
@test_approx_eq kal[:zend] alpha_hat[:, end]


# Durbin-Koopman smoother with anticipated shocks
alpha_hat, eta_hat = durbin_koopman_smoother(m, data, TTT, RRR, CCC, QQ, ZZ, DD,
    MM, EE, kal[:z0], kal[:vz0])

@test_approx_eq exp_alpha_hat alpha_hat
@test_approx_eq exp_eta_hat eta_hat

# Kalman smoother without anticipated shocks
custom_settings = Dict{Symbol, Setting}(
    :n_anticipated_shocks => Setting(:n_anticipated_shocks, 0),
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2016-Q1")))
m = Model990(custom_settings = custom_settings, testing = true)
data = data[inds_obs_no_ant(m), :]

TTT, RRR, CCC = solve(m)
meas = measurement(m, TTT, RRR, CCC)
QQ, ZZ, DD, MM, EE = meas.QQ, meas.ZZ, meas.DD, meas.MM, meas.EE
z0 = zeros(n_states_augmented(m))
P0 = QuantEcon.solve_discrete_lyapunov(TTT, RRR*QQ*RRR')

kal = kalman_filter(m, data, TTT, RRR, CCC, QQ, ZZ, DD, MM, EE, z0, P0, allout = true, include_presample = true)
alpha_hat, eta_hat = kalman_smoother(m, data, TTT, RRR, CCC, QQ, ZZ, DD, z0, P0,
    kal[:pred], kal[:vpred])
@test_matrix_approx_eq kal[:zend] alpha_hat[:, end]

# Hamilton smoother without anticipated shocks
system = compute_system(m)
alpha_hat, eta_hat = hamilton_smoother(m, data, system, kal)

@test_approx_eq kal[:zend] alpha_hat[:, end]

nothing
