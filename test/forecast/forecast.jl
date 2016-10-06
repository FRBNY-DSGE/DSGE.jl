using DSGE, HDF5, DistributedArrays
include("../util.jl")

path = dirname(@__FILE__)

# Set up arguments
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2016-Q1")),
    :use_parallel_workers => Setting(:use_parallel_workers, true),
    :n_anticipated_shocks => Setting(:n_anticipated_shocks, 6))
m = Model990(custom_settings = custom_settings)
m.testing = true

params_sim = h5open("$path/../reference/filter_args.h5","r") do h5
    read(h5, "params_sim")
end

# Add parallel workers
ndraws = 2
my_procs = addprocs(ndraws)
@everywhere using DSGE

# Set up systems
function init_systems(m, params_sim, ndraws, my_procs)
    DArray((ndraws,), my_procs, [length(my_procs)]) do I
        draw_inds = first(I)
        ndraws_local = length(draw_inds)
        localpart = Vector{System{Float64}}(ndraws_local)

        for i in draw_inds
            i_local = mod(i-1, ndraws_local) + 1

            params = squeeze(params_sim[i, :], 1)
            update!(m, params)
            localpart[i_local] = compute_system(m)
        end
        return localpart
    end
end
syses = init_systems(m, params_sim, ndraws, my_procs)

function init_states(m, syses, ndraws, my_procs)
    nstates = n_states_augmented(m)
    DArray((ndraws,), my_procs, [length(my_procs)]) do I
        [((eye(nstates) - syses[i][:TTT]) \ syses[i][:CCC])::Vector{Float64} for i in first(I)]
    end
end
z0s = init_states(m, syses, ndraws, my_procs)

# Run to compile before timing
states, observables, pseudos, shocks = DSGE.forecast(m, syses, z0s)

# Run forecast
@time states, observables, pseudos, shocks = DSGE.forecast(m, syses, z0s)

# Remove parallel workers
rmprocs(my_procs)

nothing
