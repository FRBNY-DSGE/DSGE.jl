using DSGE
using Test
using DataFrames: DataFrame
using Dates
using BenchmarkTools

@testset "Miscellaneous data handling functions" begin
    # Previous and next quarter arithmetic
    q = Dates.Date(1913,12,13)
    DSGE.prev_quarter()
    DSGE.next_quarter()
    @test DSGE.prev_quarter(q) == Dates.Date(1913,09,30)
    @test DSGE.next_quarter(q) == Dates.Date(1914,03,31)

    start_date = Dates.Date(2000,01,01)
    end_date   = Dates.Date(2010,01,01)
    DSGE.get_quarter_ends(start_date,end_date)
    DSGE.subtract_quarters(end_date, start_date)

    df1 = DataFrame(a = [1])
    global df2 = DataFrame()
    df1, df2 = DSGE.reconcile_column_names(df1, df2)
    @test names(df1) == names(df2)
    @test size(df2, 1) == 0

    quartertodate("11q3")
    quartertodate("1997q4")
    quartertodate("1985-Q1")
    @test_throws Meta.ParseError quartertodate("2005q9")
    @test_throws Meta.ParseError quartertodate("12345")

    global df = DataFrame(date = ["1913-12-23", "1992-11-14", "2002-01-01", "2014-12-19"],
                          x = [1, 2, missing, 4])
    DSGE.format_dates!(:date, df)
end

# Test hpfilter, ensuring missing data are handled correctly
missingfront = [missing; missing]
missingback  = [missing]
y  = [0.; 1.; 0.; 1.]
yn = [missingfront; y; missingback]
yf, yt   = hpfilter(y, 1600)
yfn, ytn = hpfilter(yn, 1600)

@testset "HP Filter and missing data handling" begin
    @test isequal(missingfront, yfn[1:length(missingfront)])
    @test isequal(missingfront, ytn[1:length(missingfront)])
    @test isequal(missingback, yfn[(end-length(missingback)+1):end])
    @test isequal(missingback, ytn[(end-length(missingback)+1):end])
end

################
# Benchmarking #
################
# Flip to true to run; off by default. None hit the FRED API.
run_benchmarks = false

if run_benchmarks
    bench_y     = Float64[sin(i / 4) for i in 1:280]   # ~70yr quarterly series
    bench_start = Dates.Date(2000, 01, 01)             # ~10yr span
    bench_end   = Dates.Date(2010, 01, 01)
    bench_dates_df = DataFrame(date = ["1913-12-23", "1992-11-14",
                                       "2002-01-01", "2014-12-19"])

    b_hpfilter    = @benchmark hpfilter($bench_y, 1600)
    b_quarterends = @benchmark DSGE.get_quarter_ends($bench_start, $bench_end)
    b_quartertodate = @benchmark quartertodate("1997q4")
    # evals=1: format_dates! mutates in place, so rebuild a fresh df per eval.
    b_formatdates = @benchmark DSGE.format_dates!(:date, df) setup=(df = copy($bench_dates_df)) evals=1

    println("\n===== data/misc benchmark results =====")
    for (name, b) in [("hpfilter (n=280)", b_hpfilter),
                      ("get_quarter_ends", b_quarterends),
                      ("quartertodate   ", b_quartertodate),
                      ("format_dates!   ", b_formatdates)]
        println(name, "  time:   ", BenchmarkTools.prettytime(median(b).time),
                "   memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end

nothing
