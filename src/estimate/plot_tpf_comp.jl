using ClusterManagers,HDF5,DataFrames,Plots
using QuantEcon:solve_discrete_lyapunov
using DSGE

include("plot_tpf.jl")
kalman = plot_tpf("kalman",identity=false)

procs = addprocs_sge(10,queue="background.q")
@everywhere using DSGE
@everywhere include("plot_tpf.jl")
tpfID = plot_tpf("tpf", identity=true)
rmprocs(procs,10)

procs = addprocs_sge(10,queue="background.q")
@everywhere using DSGE
@everywhere include("plot_tpf.jl")
tpfCovs = plot_tpf("tpf", identity=false)
rmprocs(procs,10)

plotly()
plot(kalman, color=:blue, linewidth=2, label="Kalman")
plot!(tpfCovs, color=:green, linewidth=2,label="cov_s")
plot!(tpfID, color=:red, linewidth=2,label="identity")
plot!(legend=:bottomright, xlabel="time", ylabel="log likelihood")
gui()