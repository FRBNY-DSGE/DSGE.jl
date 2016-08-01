# This is an example testing script to show that the code in
# forecast/filter.jl runs and produces nonempty results.  It does not
# guarantee correctness.
# TODO: compare against MATLAB for replicability

import Base.filter
using DSGE, DataFrames, HDF5
  
m = Model990()
m.testing = true
path = dirname(@__FILE__)

# Prepare data matrix
data,dates = h5open("$path/../reference/filter_args.h5","r") do h5
    read(h5,"data"), read(h5, "dates")
end

df        = DataFrame(data)  
df[:date] = Date(dates)
    
# Set up system matrices
ndraws = 2
sys = compute_system(m)
syses = repmat([sys],ndraws)

# Add parallel workers, filter, and remove parallel workers 
my_procs = addprocs(ndraws)
kals = DSGE.filter(m, df, syses, allout=true)
rmprocs(my_procs)

@time kals = DSGE.filter(m, df, syses, allout=true)

# Test that filtered states are correct


nothing