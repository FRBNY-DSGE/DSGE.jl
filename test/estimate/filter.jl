using DSGE, DataFrames, JLD2
using Dates, Test, BenchmarkTools

path = dirname(@__FILE__)
isdefined(@__MODULE__, :as_dataframe) || include(joinpath(@__DIR__, "..", "jld2_compat.jl"))

# Set up arguments
m = AnSchorfheide(testing = true)
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))

isdefined(@__MODULE__, :as_dataframe) || include(joinpath(@__DIR__, "..", "jld2_compat.jl"))

df, system, z0, P0 = JLD2.jldopen("$path/../reference/forecast_args.jld2", "r") do file
    as_dataframe(read(file, "df")), read(file, "system"), read(file, "z0"), read(file, "P0")
end
df = as_dataframe(df)

# Read expected output
exp_kal = JLD2.jldopen("$path/../reference/filter_out.jld2", "r") do file
    read(file, "exp_kal")
end
df2 = DataFrame()
df2[!, :date] = df[!, :date]
df2[!, :obs_cpi] = df[!, :obs_cpi]
df2[!, :obs_gdp] = df[!, :obs_gdp]
df2[!, :obs_nominalrate] = df[!, :obs_nominalrate]

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

################
# Benchmarking #
################
# Flip to true to run; off by default.
run_benchmarks = false

if run_benchmarks
    b_default = @benchmark DSGE.filter($m, $df, $system)
    b_init    = @benchmark DSGE.filter($m, $df, $system, $z0, $P0)

    println("\n===== estimate/filter benchmark results =====")
    for (name, b) in [("filter (default init)", b_default),
                      ("filter (z0/P0 init)  ", b_init)]
        println(name, "  time:   ", BenchmarkTools.prettytime(median(b).time),
                "   memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end

nothing
