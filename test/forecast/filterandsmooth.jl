using DSGE, DataFrames, JLD, DistributedArrays
include("../util.jl")

path = dirname(@__FILE__())

# Set up arguments
custom_settings = Dict{Symbol, Setting}(
    :n_anticipated_shocks => Setting(:n_anticipated_shocks, 6),
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")),
    :use_parallel_workers => Setting(:use_parallel_workers, true))
m = Model990(custom_settings = custom_settings, testing = true)

df, systems, z0, vz0 = jldopen("$path/../reference/forecast_args.jld","r") do file
    read(file, "df"), read(file, "systems"), read(file, "z0"), read(file, "vz0")
end

# Add parallel workers
ndraws = length(systems) # 2
my_procs = addprocs(ndraws)
@everywhere using DSGE

systems = distribute(systems; procs = my_procs, dist = [ndraws])

# Run to compile before timing
filterandsmooth_all(m, df, systems; procs = my_procs)

# Read in expected output
exp_states, exp_shocks, exp_pseudo =
    jldopen("$path/../reference/filterandsmooth_out.jld", "r") do file
        read(file, "exp_states"),
        read(file, "exp_shocks"),
        read(file, "exp_pseudo")
    end

for smoother in [:durbin_koopman, :kalman]
    m <= Setting(:forecast_smoother, smoother)

    # Without providing z0 and vz0
    @time states, shocks, pseudo = filterandsmooth_all(m, df, systems; procs = my_procs)

    @test_matrix_approx_eq exp_states[(smoother, :no_z0)] convert(Array, states)
    @test_matrix_approx_eq exp_shocks[(smoother, :no_z0)] convert(Array, shocks)
    @test_matrix_approx_eq exp_pseudo[(smoother, :no_z0)] convert(Array, pseudo)

    # Providing z0 and vz0
    @time states, shocks, pseudo = filterandsmooth_all(m, df, systems, z0, vz0; procs = my_procs)

    @test_matrix_approx_eq exp_states[(smoother, :z0)] convert(Array, states)
    @test_matrix_approx_eq exp_shocks[(smoother, :z0)] convert(Array, shocks)
    @test_matrix_approx_eq exp_pseudo[(smoother, :z0)] convert(Array, pseudo)
end

# Remove parallel workers
rmprocs(my_procs)


nothing