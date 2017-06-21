using DSGE,DataFrames,HDF5

m=AnSchorfheide(testing=true)
m<=Setting(:date_forecast_start,quartertodate("2015-Q4"))

path=dirname(@__FILE__)

df,system,s0,eps0 = jldopen($"path/..reference/forecast_args.jld","r") do file
    read(file,"df"),read(file,"system"),read(file,"s0",read(file,"eps0"
end

#read expected output
exp 