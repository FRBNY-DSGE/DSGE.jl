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

s <= Setting(:n_particles, 400)
s <= Setting(:n_Φ, 200)
s <= Setting(:λ, 2.0)
s <= Setting(:n_smc_blocks, 1)
s <= Setting(:parallel, true)
s <= Setting(:step_size_smc,0.5)
s <= Setting(:n_MH_steps_smc, 5)
s <= Setting(:resampler_smc, :multinomial)
s <= Setting(:target_accept, 0.25)
s <= Setting(:draw_type, "mode")
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

##to precompute the hessian and mode parameters
#estimate(s, data, method=:MH)
#estimate(s, data[:,100:200], method =:MH, proposal_covariance=ones(3,3))
##
smc(s,data[:,100:120],verbose=:low)

#=
@test_approx_eq resamples readcsv("resamples.csv")
@test_approx_eq wtsim readcsv("wtsim.csv")
@test_approx_eq draws readcsv("draws.csv")
=#

rm("resamples.csv")
rm("wtsim.csv")
rm("draws.csv")

nothing
