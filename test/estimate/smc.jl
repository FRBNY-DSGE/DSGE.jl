using DSGE
using HDF5, Base.Test
using DataFrames

include("../util.jl")
path = dirname(@__FILE__)

# Set up model for testing
s = AnSchorfheide()
s.testing=true
s <= Setting(:use_parallel_workers,false)

s <= Setting(:saveroot, "$path/../../save")
s <= Setting(:dataroot, "$path/../../save/input_data")
s <= Setting(:data_vintage, "160706")
s <= Setting(:date_presample_start, quartertodate("1983-Q1"))
s <= Setting(:date_mainsample_start, quartertodate("1983-Q1"))
s <= Setting(:date_mainsample_end, quartertodate("2002-Q4"))
s <= Setting(:date_zlbregime_start, quartertodate("2002-Q4"))

s <= Setting(:n_particles, 200)
s <= Setting(:n_Φ, 10)
s <= Setting(:λ, 2.0)
s <= Setting(:n_smc_blocks, 1)

try 
    rm("resamples.csv")
    rm("wtsim.csv")
    rm("draws.csv")
end

# Reading in reference 
data_mat = h5read("$path/../reference/smc.h5","data")
resamples = h5read("$path/../reference/smc.h5","resamples")
wtsim = h5read("$path/../reference/smc.h5","wtsim")
draws = h5read("$path/../reference/smc.h5","draws")

smc(s,data_mat)

@test_approx_eq resamples readcsv("resamples.csv")
@test_approx_eq wtsim readcsv("wtsim.csv")
@test_approx_eq draws readcsv("draws.csv")

rm("resamples.csv")
rm("wtsim.csv")
rm("draws.csv")

nothing
