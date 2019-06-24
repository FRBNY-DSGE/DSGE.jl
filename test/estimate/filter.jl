using DSGE, DataFrames, JLD2
using Dates, Test

path = dirname(@__FILE__)

# Set up arguments
m = AnSchorfheide(testing = true)
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))

df, system, z0, P0 = JLD2.jldopen("$path/../reference/forecast_args.jld2", "r") do file
    read(file, "df"), read(file, "system"), read(file, "z0"), read(file, "P0")
end

# Read expected output
exp_kal = JLD2.jldopen("$path/../reference/filter_out.jld2", "r") do file
    read(file, "exp_kal")
end
df2 = DataFrame()
df2[:date] = df[:date]
df2[:obs_cpi] = df[:obs_cpi]
df2[:obs_gdp] = df[:obs_gdp]
df2[:obs_nominalrate] = df[:obs_nominalrate]

# Without providing z0 and P0
@testset "Check Kalman filter outputs without initializing state/state-covariance" begin
    kal = DSGE.filter(m, df, system)
    for out in fieldnames(typeof(kal))
        global expect = exp_kal[out]
        global actual = kal[out]

        if ndims(expect) == 0
            @test expect ≈ actual
        else
            @test @test_matrix_approx_eq(expect, actual)
        end
    end
end

# Providing z0 and P0
@testset "Check Kalman filter outputs initializing state/state-covariance" begin
    kal = DSGE.filter(m, df, system, z0, P0)
    for out in fieldnames(typeof(kal))
        global expect = exp_kal[out]
        global actual = kal[out]

        if ndims(expect) == 0
            @test expect ≈ actual
        else
            @test @test_matrix_approx_eq(expect, actual)
        end
    end
end


nothing
