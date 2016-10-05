using DSGE, DataFrames, HDF5, DistributedArrays
include("../util.jl")

path = dirname(@__FILE__())

# Set up arguments
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")),
    :use_parallel_workers => Setting(:use_parallel_workers, true))
m = Model990(custom_settings = custom_settings, testing = true)

params_sim = h5open("$path/../reference/filter_args.h5","r") do h5
    read(h5, "params_sim")
end

df = load_data(m; try_disk = true, verbose = :none)

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

z0  = (eye(n_states_augmented(m)) - syses[1][:TTT]) \ syses[1][:CCC]
vz0 = QuantEcon.solve_discrete_lyapunov(syses[1][:TTT], syses[1][:RRR]*syses[1][:QQ]*syses[1][:RRR]')
kals = DSGE.filter(m, df, syses, z0, vz0; allout = true, my_procs = my_procs)

# Run to compile before timing
alpha_hats, eta_hats = smooth(m, df, syses, kals; my_procs = my_procs)

# Call smoother and test
for smoother in [:durbin_koopman, :kalman]
    m <= Setting(:forecast_smoother, smoother)

    @time alpha_hats, eta_hats = smooth(m, df, syses, kals; my_procs = my_procs)

    exp_alpha_hats = Vector{Matrix{Float64}}(ndraws)
    exp_eta_hats   = Vector{Matrix{Float64}}(ndraws)
    for i = 1:ndraws
        exp_alpha_hats[i], exp_eta_hats[i] = if forecast_smoother(m) == :durbin_koopman
            durbin_koopman_smoother(m, df, syses[i], kals[i][:z0], kals[i][:vz0])
        elseif forecast_smoother(m) == :kalman
            kalman_smoother(m, df, syses[i], kals[i][:z0], kals[i][:vz0], kals[i][:pred], kals[i][:vpred])
        end

        @test_matrix_approx_eq exp_alpha_hats[i] convert(Array, slice(alpha_hats, i, :, :))
        @test_matrix_approx_eq exp_eta_hats[i]   convert(Array, slice(eta_hats, i, :, :))
    end
end

# Remove parallel workers
rmprocs(my_procs)

nothing