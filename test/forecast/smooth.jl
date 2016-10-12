using DSGE, DataFrames, HDF5, DistributedArrays
include("../util.jl")

path = dirname(@__FILE__())

# Set up arguments
custom_settings = Dict{Symbol, Setting}(
    :n_anticipated_shocks => Setting(:n_anticipated_shocks, 6),
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")),
    :use_parallel_workers => Setting(:use_parallel_workers, true))
m = Model990(custom_settings = custom_settings, testing = true)

df, systems, kals = jldopen("$path/../reference/forecast_args.jld","r") do file
    read(file, "df"), read(file, "systems"), read(file, "kals")
end

# Add parallel workers
ndraws = length(systems) # 2
my_procs = addprocs(ndraws)
@everywhere using DSGE

systems = distribute(systems; procs = my_procs, dist = [ndraws])
kals    = distribute(kals;    procs = my_procs, dist = [ndraws])

# Run to compile before timing
states, shocks = smooth_all(m, df, systems, kals; procs = my_procs)

# Read expected output
exp_states, exp_shocks = jldopen("$path/../reference/smooth_out.jld", "r") do file
    read(file, "exp_states"), read(file, "exp_shocks")
end

# Call smoother and test
for smoother in [:durbin_koopman, :kalman]
    m <= Setting(:forecast_smoother, smoother)

    @time states, shocks = smooth_all(m, df, systems, kals; procs = my_procs)

    @test_matrix_approx_eq exp_states[smoother] convert(Array, states)
    @test_matrix_approx_eq exp_shocks[smoother] convert(Array, shocks)
end

# Remove parallel workers
rmprocs(my_procs)

nothing