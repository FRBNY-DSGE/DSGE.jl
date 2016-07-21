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

alpha_hat, eta_hat = kalman_smoother(A0, P0, y, pred, vpred, T, R, Q, Z, b,
                                     nant, antlags, peachcount, psize)

exp_alpha_hat, exp_eta_hat =
    h5open("$path/../reference/kalman_smoother_out.h5", "r") do h5
        read(h5, "alpha_hat"), read(h5, "eta_hat")
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

TTT, RRR, CCC = solve(m)
meas = measurement(m, TTT, RRR, CCC)
QQ, ZZ, DD = meas.QQ, meas.ZZ, meas.DD

alpha_hat, eta_hat = durbin_koopman_smoother(m, data, TTT, RRR, CCC, QQ, ZZ, DD, P0)

# exp_alpha_hat, exp_eta_hat =
#     h5open("$path/../reference/durbin_koopman_smoother_out.h5", "r") do h5
#         read(h5, "alpha_hat"), read(h5, "eta_hat")
# end

# @test_approx_eq exp_alpha_hat alpha_hat
# @test_approx_eq exp_eta_hat eta_hat


nothing
