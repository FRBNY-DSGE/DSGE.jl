using DSGE
using HDF5

path = dirname(@__FILE__)

# 
m = Model990()
m.testing = true

df = load_data(m; try_disk=true, verbose=:none)
forecast_one(m, df; input_type=:init, output_type=:forecast, cond_type=:none)

nothing
