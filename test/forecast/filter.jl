using DSGE, DataFrames, JLD
include("../util.jl")

path = dirname(@__FILE__)

# Set up arguments
custom_settings = Dict{Symbol, Setting}(
    :n_anticipated_shocks => Setting(:n_anticipated_shocks, 6),
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")),
    :use_parallel_workers => Setting(:use_parallel_workers, true))
m = Model990(custom_settings = custom_settings, testing = true)

df, system, z0, vz0 = jldopen("$path/../reference/forecast_args.jld", "r") do file
    read(file, "df"), read(file, "system"), read(file, "z0"), read(file, "vz0")
end

# Read expected output
exp_kals_no_z0, exp_kals_z0 = jldopen("$path/../reference/filter_out.jld","r") do file
    read(file, "exp_kals_no_z0"), read(file, "exp_kals_z0")
end

# Without providing z0 and vz0
kal = DSGE.filter(m, df, system; allout = true)
for out in fieldnames(kal)
    expect = exp_kals_no_z0[1][out]
    actual = kal[out]

    ndims(expect) == 0 ? @test_approx_eq(expect, actual) : @test_matrix_approx_eq(expect, actual)
end

# Providing z0 and vz0
kal = DSGE.filter(m, df, system, z0, vz0; allout = true)
for out in fieldnames(kal)
    expect = exp_kals_z0[1][out]
    actual = kal[out]

    ndims(expect) == 0 ? @test_approx_eq(expect, actual) : @test_matrix_approx_eq(expect, actual)
end


nothing