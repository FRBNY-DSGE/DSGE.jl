using BenchmarkTools
using ClusterManagers
using HDF5, JLD, Base.Test
using DataFrames
using DSGE

m = AnSchorfheide()

data = h5read("smc.h5","data")

