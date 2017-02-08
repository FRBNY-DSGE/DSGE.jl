using DSGE, DataFrames, JLD
include("../util.jl")

path = dirname(@__FILE__)

# Set up arguments
m = AnSchorfheide(testing = true)
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))

df, system, z0, vz0 = jldopen("$path/../reference/forecast_args.jld", "r") do file
    read(file, "df"), read(file, "system"), read(file, "z0"), read(file, "vz0")
end

# Read expected output
exp_kal_no_z0, exp_kal_z0 = jldopen("$path/../reference/filter_out.jld", "r") do file
    read(file, "exp_kal_no_z0"), read(file, "exp_kal_z0")
end

# Without providing z0 and vz0
kal = DSGE.filter(m, df, system; allout = true)
for out in fieldnames(kal)
    expect = exp_kal_no_z0[out]
    actual = kal[out]

    ndims(expect) == 0 ? @test_approx_eq(expect, actual) : @test_matrix_approx_eq(expect, actual)
end

# Providing z0 and vz0
kal = DSGE.filter(m, df, system, z0, vz0; allout = true)
for out in fieldnames(kal)
    expect = exp_kal_z0[out]
    actual = kal[out]

    ndims(expect) == 0 ? @test_approx_eq(expect, actual) : @test_matrix_approx_eq(expect, actual)
end

nothing