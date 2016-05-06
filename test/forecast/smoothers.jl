using HDF5, Base.Test
include("../util.jl")

path = dirname(@__FILE__)

# Initialize arguments to function
h5 = h5open("$path/../reference/kalman_smoother_args.h5")
for arg in ["A0", "C", "P0", "Q", "R", "T", "Z", "antlags", "b", "nant", "peachcount", "pred", "psize", "vpred", "y"]
    eval(parse("$arg = read(h5, \"$arg\")"))
end
for arg in ["b"]
    eval(parse("$arg = reshape(read(h5, \"$arg\"), length($arg), 1)"))
end
    
for arg in ["nant", "antlags","psize", "peachcount"]
    eval(parse("$arg = round(Int, $arg[1])"))
end


# Method with all arguments provided (9)
kalsmooth = kalman_smoother(A0, P0, y, pred, vpred, T, R, Q, Z, b, nant, antlags, peachcount, psize)
alpha_hat = kalsmooth.states
eta_hat   = kalsmooth.shocks 

exp_alpha_hat, exp_eta_hat = h5open("$path/../reference/kalman_smoother_out.h5", "r") do h5
    read(h5,"alpha_hat"), read(h5,"eta_hat")
end

@test_approx_eq exp_alpha_hat alpha_hat
@test_approx_eq exp_eta_hat eta_hat


close(h5)


nothing
