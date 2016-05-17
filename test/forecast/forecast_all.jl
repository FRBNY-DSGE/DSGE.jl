using DSGE
using HDF5

path = dirname(@__FILE__)

# 
m = Model990()

data = load_data(m)
forecast_one(m, data; input_type=:mode, output_type=:forecast, cond_type=:none)

nothing
