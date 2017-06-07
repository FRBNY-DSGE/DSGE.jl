using DSGE
using HDF5, Base.Test
using DataFrames

path = dirname(@__FILE__)

# Set up model for testing

custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")))
s = AnSchorfheide(custom_settings = custom_settings, testing = true)

s <= Setting(:saveroot, "$path/../../save")
s <= Setting(:dataroot, "$path/../../save/input_data")

s <= Setting(:n_particles, 200)
s <= Setting(:n_Φ, 10)
s <= Setting(:λ, 1.)
s <= Setting(:n_smc_blocks, 1)

try
    rm("resamples.csv")
    rm("wtsim.csv")
    rm("draws.csv")
end

# Reading in reference
data = h5read("$path/../reference/smc.h5","data")
resamples = h5read("$path/../reference/smc.h5","resamples")
wtsim = h5read("$path/../reference/smc.h5","wtsim")
draws = h5read("$path/../reference/smc.h5","draws")

smc(s,data)

@test_approx_eq resamples readcsv("resamples.csv")
@test_approx_eq wtsim readcsv("wtsim.csv")
@test_approx_eq draws readcsv("draws.csv")

rm("resamples.csv")
rm("wtsim.csv")
rm("draws.csv")

nothing
