using DSGE, HDF5
include("../util.jl")

path = dirname(@__FILE__)

# Set up arguments
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2016-Q1")),
    :use_parallel_workers => Setting(:use_parallel_workers, true),
    :n_anticipated_shocks => Setting(:n_anticipated_shocks, 6))
m = Model990(custom_settings = custom_settings)
m.testing = true

eta_hat = h5open("$path/../reference/durbin_koopman_smoother_out.h5", "r") do h5
    read(h5, "eta_hat")
end

ndraws = 2
syses = Vector{System{Float64}}(ndraws)
histshocks = Vector{Matrix{Float64}}(ndraws)
for i = 1:ndraws
    syses[i] = compute_system(m)
    histshocks[i] = eta_hat
end

# Add parallel workers
my_procs = addprocs(ndraws)
@everywhere using DSGE
states, observables, pseudos = DSGE.shock_decompositions(m, syses, histshocks)

# Run forecast
@time states, observables, pseudos = DSGE.shock_decompositions(m, syses, histshocks)

# Remove parallel workers
rmprocs(my_procs)

nothing
