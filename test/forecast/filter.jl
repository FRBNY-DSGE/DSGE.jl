import Base.filter
using DSGE, DataFrames, HDF5
include("../util.jl")

path = dirname(@__FILE__)

# Set up arguments
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")),
    :use_parallel_workers => Setting(:use_parallel_workers, true))
m = Model990(custom_settings = custom_settings, testing = true)

params_sim = h5open("$path/../reference/filter_args.h5","r") do h5
    read(h5, "params_sim")
end

df = load_data(m; try_disk = true, verbose = :none)

ndraws = 2
syses = Vector{System{Float64}}(ndraws)
for i = 1:ndraws
    params = squeeze(params_sim[i, :], 1)
    update!(m, params)
    syses[i] = compute_system(m)
end

z0  = (eye(n_states_augmented(m)) - syses[1][:TTT]) \ syses[1][:CCC]
vz0 = QuantEcon.solve_discrete_lyapunov(syses[1][:TTT], syses[1][:RRR]*syses[1][:QQ]*syses[1][:RRR]')

# Add parallel workers
my_procs = addprocs(ndraws)
@everywhere using DSGE
kals = DSGE.filter(m, df, syses; allout = true)
kals = DSGE.filter(m, df, syses, z0, vz0; allout = true)

# Without providing z0 and vz0
@time kals = DSGE.filter(m, df, syses; allout = true)

exp_kals = Vector{DSGE.Kalman{Float64}}(ndraws)
for i = 1:ndraws
    exp_kals[i] = kalman_filter(m, df_to_matrix(m, df), syses[i][:TTT], syses[i][:CCC], syses[i][:ZZ], syses[i][:DD], syses[i][:VVall]; allout = true)
end

for i = 1:ndraws
    for out in fieldnames(kals[1])
        expect = exp_kals[i][out]
        actual = kals[i][out]

        if ndims(expect) == 0
            @test_approx_eq expect actual
        else
            @test_matrix_approx_eq expect actual
        end
    end
end

# Providing z0 and vz0
@time kals = DSGE.filter(m, df, syses, z0, vz0; allout = true)

exp_kals = Vector{DSGE.Kalman{Float64}}(ndraws)
for i = 1:ndraws
    exp_kals[i] = kalman_filter(m, df_to_matrix(m, df), syses[i][:TTT], syses[i][:CCC], syses[i][:ZZ], syses[i][:DD], syses[i][:VVall], z0, vz0; allout = true)
end

for i = 1:ndraws
    for out in fieldnames(kals[1])
        expect = exp_kals[i][out]
        actual = kals[i][out]

        if ndims(expect) == 0
            @test_approx_eq expect actual
        else
            @test_matrix_approx_eq expect actual
        end
    end
end

# Remove parallel workers
rmprocs(my_procs)

nothing