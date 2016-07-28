using DSGE, HDF5

path = dirname(@__FILE__())

m = Model990()
data = h5open("$path/../reference/filter_args.h5","r") do h5
    read(h5,"data")
end

sys = compute_system(m)

alpha_hat, eta_hat = DSGE.filterandsmooth(m, data, sys, allout=true)

nothing