using DSGE
using HDF5, Base.Test
using DataFrames

path = dirname(@__FILE__)

# Set up model for testing

custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")))
s = AnSchorfheide(custom_settings = custom_settings, testing = true)

saveroot = normpath(joinpath(dirname(@__FILE__),"..","..","save"))
dataroot = normpath(joinpath(dirname(@__FILE__),"..","..","save","input_data"))

s <= Setting(:saveroot, saveroot)
s <= Setting(:dataroot, dataroot)
s <= Setting(:n_particles, 400)
s <= Setting(:n_Φ, 200)
s <= Setting(:λ, 2.0)
s <= Setting(:n_smc_blocks, 1)
s <= Setting(:use_parallel_workers, false)
s <= Setting(:step_size_smc,0.5)
s <= Setting(:n_MH_steps_smc, 5)
s <= Setting(:resampler_smc, :systematic)
s <= Setting(:target_accept, 0.25)

data = h5read("$path/../reference/smc.h5","data")
cloud = ParticleCloud(s,get_setting(s,:n_particles))

# Testing initial_draw()
s <= Setting(:initial_draw_source, :normal)
initial_draw(s, data, cloud)
@assert sum(!isfinite(DSGE.get_logpost(cloud))) == 0
@assert sum(!isfinite(DSGE.get_loglh(cloud))) == 0

s <= Setting(:initial_draw_source, :prior)
initial_draw(s, data, cloud)
@assert sum(!isfinite(DSGE.get_logpost(cloud))) == 0
@assert sum(!isfinite(DSGE.get_loglh(cloud))) == 0

# Reading in reference
data = h5read("$path/../reference/smc.h5","data")
weights = h5read("$path/../reference/smc.h5","weights")
draws = h5read("$path/../reference/smc.h5","draws")
mutated_values = h5read("$path/../reference/smc.h5","mutated_values")

# to precompute the hessian and mode parameters
smc(s,data,verbose=:none)

@test_approx_eq draws readcsv("draws.csv")
@test_approx_eq weights readcsv("weights.csv")
@test_approx_eq mutated_values readcsv("mutated_values.csv")

rm("draws.csv")
rm("weights.csv")
rm("mutated_values.csv")

nothing
