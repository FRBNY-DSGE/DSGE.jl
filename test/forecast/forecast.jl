using DSGE
include("../util.jl")

path = dirname(@__FILE__)

# Set up arguments
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2016-Q1")),
    :use_parallel_workers => Setting(:use_parallel_workers, true))
m = Model990(custom_settings = custom_settings)
m.testing = true

params_sim = h5open("$path/../reference/filter_args.h5","r") do h5
    read(h5, "params_sim")
end

ndraws = 2
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
