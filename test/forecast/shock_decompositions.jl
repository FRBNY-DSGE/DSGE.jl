using DSGE, HDF5, DistributedArrays
include("../util.jl")

path = dirname(@__FILE__)

# Set up arguments
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2016-Q1")),
    :date_forecast_end    => Setting(:date_forecast_end, quartertodate("2016-Q2")),
    :use_parallel_workers => Setting(:use_parallel_workers, true),
    :n_anticipated_shocks => Setting(:n_anticipated_shocks, 6),
    :forecast_pseudoobservables => Setting(:forecast_pseudoobservables, true))
m = Model990(custom_settings = custom_settings, testing = true)

eta_hat = h5open("$path/../reference/durbin_koopman_smoother_out.h5", "r") do h5
    read(h5, "eta_hat")
end

# Add parallel workers
ndraws = 2
my_procs = addprocs(ndraws)
@everywhere using DSGE

system = compute_system(m)
systems = dfill(system, (ndraws,), my_procs, [ndraws])
histshocks = repeat(reshape(eta_hat, (1, size(eta_hat)...)), outer = [ndraws, 1, 1])
histshocks = distribute(histshocks; procs = my_procs, dist = [ndraws, 1, 1])

# Run to compile before timing
states, observables, pseudos = DSGE.shock_decompositions(m, systems, histshocks)

# Run shock decompositions
@time states, observables, pseudos = DSGE.shock_decompositions(m, systems, histshocks)

@assert !isnull(shockdec_startdate(m))
nperiods = DSGE.subtract_quarters(date_forecast_end(m), get(shockdec_startdate(m))) + 1
nstates  = n_states_augmented(m)
nobs     = n_observables(m)
npseudo  = n_pseudoobservables(m)
nshocks  = n_shocks_exogenous(m)

for i = 1:ndraws
    @assert size(slice(states,      i, :, :, :)) == (nstates, nperiods, nshocks)
    @assert size(slice(observables, i, :, :, :)) == (nobs,    nperiods, nshocks)
    @assert size(slice(pseudos,     i, :, :, :)) == (npseudo, nperiods, nshocks)
end

# Run forecast again, with shockdec_startdate null
m <= Setting(:shockdec_startdate, Nullable{Date}())
@time states, observables, pseudos = DSGE.shock_decompositions(m, systems, histshocks)

nperiods = DSGE.subtract_quarters(date_forecast_end(m), date_prezlb_start(m)) + 1
for i = 1:ndraws
    @assert size(slice(states,      i, :, :, :)) == (nstates, nperiods, nshocks)
    @assert size(slice(observables, i, :, :, :)) == (nobs,    nperiods, nshocks)
    @assert size(slice(pseudos,     i, :, :, :)) == (npseudo, nperiods, nshocks)
end

# Remove parallel workers
rmprocs(my_procs)

nothing
