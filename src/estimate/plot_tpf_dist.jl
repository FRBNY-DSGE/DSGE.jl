using ClusterManagers,Plots

path = dirname(@__FILE__)

addprocs_sge(10, queue="background.q")
@everywhere using DSGE, HDF5, DataFrames
@everywhere using QuantEcon: solve_discrete_lyapunov
@everywhere include("tpf_error.jl")
error = zeros(5)

error = pmap(i -> tpf_error(), 1:5)

plotly()
histogram(error)
gui()

# The following code is still subject to the increasing runtime issue despite removing and re-adding all workers
# for i=1:30
#     @show i
#     procs = addprocs_sge(10,queue="background.q")
#     @everywhere using DSGE
#     @everywhere include("tpf_error.jl")
#     error[i] = tpf_error()
#     h5open("$path/../../test/reference/error.h5","w") do file
#         write(file,"error",error)
#     end
#     rmprocs(procs,10)
# end