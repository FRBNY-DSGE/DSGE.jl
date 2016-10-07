using DSGE, HDF5, JLD, DistributedArrays

ndraws = 2

# Set up arguments
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")),
    :use_parallel_workers => Setting(:use_parallel_workers, true),
    :n_anticipated_shocks => Setting(:n_anticipated_shocks, 6),
    :forecast_pseudoobservables => Setting(:forecast_pseudoobservables, true))
m = Model990(custom_settings = custom_settings, testing = true)
ndraws = 2

df = load_data(m; try_disk = true, verbose = :none)

# Set up systems and states
systems = Vector{System{Float64}}(ndraws)
states = Vector{Vector{Float64}}(ndraws)

params_sim = h5open("$path/../reference/filter_args.h5","r") do h5
    read(h5, "params_sim")
end
nstates = n_states_augmented(m)

for i = 1:ndraws
    params = squeeze(params_sim[i, :], 1)
    update!(m, params)
    systems[i] = compute_system(m)
    states[i]  = zeros(nstates)
end

# Initial state and covariance
z0  = zeros(nstates)
vz0 = QuantEcon.solve_discrete_lyapunov(systems[1][:TTT], systems[1][:RRR]*systems[1][:QQ]*systems[1][:RRR]')
ind_ant1 = m.endogenous_states[:rm_tl1]
ind_antn = m.endogenous_states[symbol("rm_tl$(n_anticipated_shocks(m))")]
shock_inds_ant = ind_ant1:ind_antn
vz0[:, shock_inds_ant] = 0;
vz0[shock_inds_ant, :] = 0;

# Add parallel workers
ndraws = length(systems) # 2
my_procs = addprocs(ndraws)
@everywhere using DSGE

# Kalman objects
systems_dist = distribute(systems; procs = my_procs, dist = [ndraws])
kals = DSGE.filter(m, df, systems_dist, z0, vz0; allout = true, procs = my_procs)
kals = convert(Array, kals)

# Smoothed shocks
kals_dist = distribute(kals; procs = my_procs, dist = [ndraws])
_, histshocks = smooth(m, df, systems_dist, kals_dist; procs = my_procs)
histshocks = convert(Array, histshocks)

# Remove parallel workers
rmprocs(my_procs)

# Save to JLD
jldopen("$path/../reference/forecast_args.jld", "w") do file
    write(file, "df", df),
    write(file, "systems", systems),
    write(file, "states", states),
    write(file, "z0", z0),
    write(file, "vz0", vz0),
    write(file, "kals", kals),
    write(file, "histshocks", histshocks)
end

nothing