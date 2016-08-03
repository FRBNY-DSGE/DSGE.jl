using DSGE, DataFrames, HDF5
include("../util.jl")

path = dirname(@__FILE__())

# Set up arguments
m = Model990()
m.testing = true
m <= Setting(:date_forecast_start, quartertodate("2016-Q1"))
m <= Setting(:use_parallel_workers, true)

data, dates, params_sim = h5open("$path/../reference/filter_args.h5","r") do h5
    read(h5, "data"), read(h5, "dates"), read(h5, "params_sim")
end

df        = DataFrame(data)
df[:date] = Date(dates)

ndraws = 2
syses = Vector{System{Float64}}(ndraws)
for i = 1:ndraws
    params = squeeze(params_sim[i, :], 1)
    update!(m, params)
    syses[i] = compute_system(m)
end

z0  = (eye(n_states_augmented(m)) - syses[1][:TTT]) \ syses[1][:CCC]
vz0 = QuantEcon.solve_discrete_lyapunov(syses[1][:TTT], syses[1][:RRR]*syses[1][:QQ]*syses[1][:RRR]')
kals = DSGE.filter(m, df, syses, z0, vz0; allout = true)

# Add parallel workers
my_procs = addprocs(ndraws)
@everywhere using DSGE
alpha_hats, eta_hats = DSGE.smooth(m, df, syses, kals)

# Call smoother and test
for smoother in [:durbin_koopman, :kalman]
    m <= Setting(:forecast_smoother, smoother)

    @time alpha_hats, eta_hats = DSGE.smooth(m, df, syses, kals)

    exp_alpha_hats = Vector{Matrix{Float64}}(ndraws)
    exp_eta_hats   = Vector{Matrix{Float64}}(ndraws)
    for i = 1:ndraws
        exp_alpha_hats[i], exp_eta_hats[i] = if forecast_smoother(m) == :durbin_koopman
            durbin_koopman_smoother(m, df, syses[i], kals[i][:z0], kals[i][:vz0])
        elseif forecast_smoother(m) == :kalman
            kalman_smoother(m, df, syses[i], kals[i][:z0], kals[i][:vz0], kals[i][:pred], kals[i][:vpred])
        end

        @test_matrix_approx_eq exp_alpha_hats[i] alpha_hats[i]
        @test_matrix_approx_eq exp_eta_hats[i] eta_hats[i]
    end
end

# Remove parallel workers
rmprocs(my_procs)

nothing