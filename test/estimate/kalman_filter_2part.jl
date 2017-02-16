using DSGE, HDF5, JLD, Base.Test
include("../util.jl")

path = dirname(@__FILE__)

# Set up
data, z0, P0 =
    h5open("$path/../reference/kalman_filter_2part_args.h5", "r") do file
        read(file, "data"), read(file, "A0"), read(file, "P0")
    end

custom_settings = Dict{Symbol, Setting}(
    :n_anticipated_shocks => Setting(:n_anticipated_shocks, 6),
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2016-Q1")))
m = Model990(custom_settings = custom_settings, testing = true)

# Kalman filter with z0 and P0 provided
kal1 = kalman_filter(m, data, z0, P0; allout = true, include_presample = true)

# Kalman filter without z0 and P0
kal2 = kalman_filter(m, data; allout = true, include_presample = true)

# Test against expected output
exp_kal = jldopen("$path/../reference/kalman_filter_2part_out.jld", "r") do file
    read(file, "exp_kal")
end

for out in [:L, :zend, :Pend, :pred, :vpred, :yprederror, # :ystdprederror, :rmsd
            :rmse, :filt, :vfilt, :z0, :vz0]
    @test_approx_eq exp_kal[out] kal1[out]
    @test_approx_eq exp_kal[out] kal2[out]
end

nothing