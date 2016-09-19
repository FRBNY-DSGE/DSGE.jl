using DSGE, DataFrames, HDF5
include("../util.jl")

path = dirname(@__FILE__())

# Set up arguments
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2016-Q1")),
    :use_parallel_workers => Setting(:use_parallel_workers, true))
m = Model990(custom_settings = custom_settings, testing = true)

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

# Add parallel workers
my_procs = addprocs(ndraws)
@everywhere using DSGE
states, shocks, pseudo = filterandsmooth(m, df, syses; allout = true)
states, shocks, pseudo = filterandsmooth(m, df, syses, z0, vz0; allout = true)

for smoother in [:durbin_koopman, :kalman]
    m <= Setting(:forecast_smoother, smoother)

    # Without providing z0 and vz0
    @time states, shocks, pseudo = filterandsmooth(m, df, syses; allout = true)

    exp_states = Vector{Matrix{Float64}}(ndraws)
    exp_shocks = Vector{Matrix{Float64}}(ndraws)
    for i = 1:ndraws
        kal = kalman_filter(m, df_to_matrix(m, df), syses[i][:TTT], syses[i][:CCC], syses[i][:ZZ],
                            syses[i][:DD], syses[i][:VVall]; allout = true)

        exp_states[i], exp_shocks[i] = if forecast_smoother(m) == :durbin_koopman
            durbin_koopman_smoother(m, df, syses[i], kal[:z0], kal[:vz0])
        elseif forecast_smoother(m) == :kalman
            kalman_smoother(m, df, syses[i], kal[:z0], kal[:vz0], kal[:pred], kal[:vpred])
        end

        @test_matrix_approx_eq exp_states[i] states[:, :, i]
        @test_matrix_approx_eq exp_shocks[i] shocks[:, :, i]
        @test all(x -> x == 0, pseudo[:, :, i])
    end

    # Providing z0 and vz0
    @time states, shocks, pseudo = filterandsmooth(m, df, syses, z0, vz0; allout = true)

    exp_states = Vector{Matrix{Float64}}(ndraws)
    exp_shocks = Vector{Matrix{Float64}}(ndraws)
    for i = 1:ndraws
        kal = kalman_filter(m, df_to_matrix(m, df), syses[i][:TTT], syses[i][:CCC], syses[i][:ZZ],
                            syses[i][:DD], syses[i][:VVall], z0, vz0; allout = true)

        exp_states[i], exp_shocks[i] = if forecast_smoother(m) == :durbin_koopman
            durbin_koopman_smoother(m, df, syses[i], kal[:z0], kal[:vz0])
        elseif forecast_smoother(m) == :kalman
            kalman_smoother(m, df, syses[i], kal[:z0], kal[:vz0], kal[:pred], kal[:vpred])
        end

        @test_matrix_approx_eq exp_states[i] states[:, :, i]
        @test_matrix_approx_eq exp_shocks[i] shocks[:, :, i]
        @test all(x -> x == 0, pseudo[:, :, i])
    end

end

# Remove parallel workers
rmprocs(my_procs)

nothing