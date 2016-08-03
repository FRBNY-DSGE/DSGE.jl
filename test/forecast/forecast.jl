using DSGE
include("../util.jl")

path = dirname(@__FILE__)

# Set up arguments
m = Model990()
m.testing = true
m <= Setting(:date_forecast_start, quartertodate("2016-Q1"))
m <= Setting(:use_parallel_workers, true)

params_sim = h5open("$path/../reference/filter_args.h5","r") do h5
    read(h5, "params_sim")
end

ndraws = size(params_sim, 1)
syses = Vector{System{Float64}}(ndraws)
z0s = Vector{Vector{Float64}}(ndraws)
for i = 1:ndraws
    params = squeeze(params_sim[i, :], 1)
    update!(m, params)
    syses[i] = compute_system(m)
    z0s[i] = (eye(n_states_augmented(m)) - syses[i][:TTT]) \ syses[i][:CCC]
end

# Add parallel workers
my_procs = addprocs(ndraws)
@everywhere using DSGE

# Run forecast
states, observables, pseudos = DSGE.forecast(m, syses, z0s)

# Remove parallel workers
rmprocs(my_procs)

nothing
