using DSGE
using Base.Test, DataFrames, HDF5, JLD
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
    :n_anticipated_shocks       => Setting(:n_anticipated_shocks, 6),
    :date_forecast_start        => Setting(:date_forecast_start, quartertodate("2016-Q1")),
    :forecast_pseudoobservables => Setting(:forecast_pseudoobservables, false))
m = Model990(custom_settings = custom_settings, testing = true)

trans  = Transition(TTT, RRR, CCC)
meas   = measurement(m, TTT, RRR, CCC)
system = System(trans, meas)

# Create DataFrame out of data matrix
obs_names = sort(collect(keys(m.observables)), by = x -> m.observables[x])
df = DataFrame()
for (var, ind) in zip(obs_names, 1:n_observables(m))
    df[var] = squeeze(data[ind, :], 1)
end
df[:date] = DSGE.quarter_range(date_presample_start(m), date_mainsample_end(m))

# Read expected output
exp_alpha_hat, exp_eta_hat = h5open("$path/../reference/kalman_smoother_out.h5", "r") do h5
    read(h5, "alpha_hat"), read(h5, "eta_hat")
end


# All smoothers with anticipated shocks
for smoother in [:hamilton, :koopman, :carter_kohn, :durbin_koopman]
    m <= Setting(:forecast_smoother, smoother)
    alpha_hat, eta_hat = smooth(m, df, system, kal)

    if smoother in [:koopman, :durbin_koopman]
        @test_approx_eq exp_alpha_hat alpha_hat
        @test_approx_eq exp_eta_hat   eta_hat
    else
        @test_approx_eq_eps exp_alpha_hat alpha_hat 2e-2
        @test_approx_eq_eps exp_eta_hat   eta_hat   1e-3
    end

    @test_approx_eq kal[:zend]    alpha_hat[:, end]
end


# Koopman smoother without anticipated shocks
custom_settings = Dict{Symbol, Setting}(
    :n_anticipated_shocks => Setting(:n_anticipated_shocks, 0),
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2016-Q1")))
m = Model990(custom_settings = custom_settings, testing = true)
data = data[inds_obs_no_ant(m), :]

system = compute_system(m)
z0 = zeros(n_states_augmented(m))
P0 = QuantEcon.solve_discrete_lyapunov(system[:TTT], system[:RRR]*system[:QQ]*system[:RRR]')

m <= Setting(:forecast_smoother, :koopman)
kal = DSGE.filter(m, data, system, z0, P0, include_presample = true)
alpha_hat, eta_hat = smooth(m, df, system, kal)
@test_matrix_approx_eq kal[:zend] alpha_hat[:, end]

nothing
