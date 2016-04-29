using DSGE
using HDF5

path = dirname(@__FILE__)

# 
m = Model990()
#m.testing = true

#datapath = inpath(m, "data", "data_REF.h5")
datapath = abspath("$path/../reference/data/data_REF.h5")
data = h5open(datapath) do file
    read(file, "data")
end

forecast_one(m, data; input_type=:mode, output_type=:forecast, cond_type=:none)

nothing
