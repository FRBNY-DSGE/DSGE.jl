using DSGE, DataFrames, JLD
using Base.Test

path = dirname(@__FILE__)

# Set up arguments
m = AnSchorfheide(testing = true)
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))

df, system, s_0, P_0 = jldopen("$path/../reference/forecast_args.jld", "r") do file
    read(file, "df"), read(file, "system"), read(file, "z0"), read(file, "P0")
end

# Read expected output
exp_kal = jldopen("$path/../reference/filter_out.jld", "r") do file
    read(file, "exp_kal")
end

# Without providing s_0 and P_0
kal = DSGE.filter(m, df, system)
for out in fieldnames(kal)
    expect = exp_kal[out]
    actual = kal[out]

    if ndims(expect) == 0
        @test expect ≈ actual
    else
        @test_matrix_approx_eq(expect, actual)
    end
end

# Providing s_0 and P_0
kal = DSGE.filter(m, df, system, s_0, P_0)
for out in fieldnames(kal)
    expect = exp_kal[out]
    actual = kal[out]

    if ndims(expect) == 0
        @test expect ≈ actual
    else
        @test_matrix_approx_eq(expect, actual)
    end
end


nothing
