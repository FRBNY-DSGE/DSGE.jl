import Base.filter
using DSGE, DataFrames, JLD, DistributedArrays
include("../util.jl")

path = dirname(@__FILE__)

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
DSGE.filter_all(m, df, systems; allout = true, procs = my_procs)
DSGE.filter_all(m, df, systems, z0, vz0; allout = true, procs = my_procs)

# Read expected output
exp_kals_no_z0, exp_kals_z0 = jldopen("$path/../reference/filter_out.jld","r") do file
    read(file, "exp_kals_no_z0"), read(file, "exp_kals_z0")
end

# Without providing z0 and vz0
@time kals = DSGE.filter_all(m, df, systems; allout = true, procs = my_procs)
for i = 1:ndraws
    for out in fieldnames(kals[1])
        expect = exp_kals_no_z0[i][out]
        actual = kals[i][out]

        ndims(expect) == 0 ? @test_approx_eq(expect, actual) : @test_matrix_approx_eq(expect, actual)
    end
end

# Providing z0 and vz0
@time kals = DSGE.filter_all(m, df, systems, z0, vz0; allout = true, procs = my_procs)
for i = 1:ndraws
    for out in fieldnames(kals[1])
        expect = exp_kals_z0[i][out]
        actual = kals[i][out]

        ndims(expect) == 0 ? @test_approx_eq(expect, actual) : @test_matrix_approx_eq(expect, actual)
    end
end

# Remove parallel workers
rmprocs(my_procs)

nothing