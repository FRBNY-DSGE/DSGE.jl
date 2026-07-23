using DSGE, Dates, BenchmarkTools, DataFrames, OrderedCollections
path = dirname(@__FILE__)

mb_empty = MeansBands()
@show mb_empty
@test isempty(mb_empty)

# Build mb_full in-memory: 3 observables, 60 quarters from 2015-Q4 to 2030-Q3,
# 5 density bands (50%–90%), matching what the stale MeansBands.jld2 reference contained.
let
    n    = 60
    vars = [:obs_gdp, :obs_cpi, :obs_nominalrate]
    start_date = quartertodate("2015-Q4")   # Date(2015, 12, 31)
    dates = [start_date + Dates.Month(3*(i-1)) for i in 1:n]
    # dates[1]   == quartertodate("2015-Q4") == Date(2015, 12, 31)
    # dates[end] == quartertodate("2030-Q3") == Date(2030, 9, 30)

    band_pcts = ["90.0%", "80.0%", "70.0%", "60.0%", "50.0%"]
    band_cols = vcat([Symbol(p * " LB") for p in band_pcts],
                     [Symbol(p * " UB") for p in band_pcts])

    means_df   = DataFrame(:date => dates, (v => zeros(n) for v in vars)...)
    bands_dict = Dict{Symbol,DataFrame}(
        v => DataFrame(:date => dates, (c => zeros(n) for c in band_cols)...)
        for v in vars
    )
    metadata = Dict{Symbol,Any}(
        :product         => :forecast,
        :class           => :obs,
        :input_type      => :full,
        :cond_type       => :none,
        :para            => :full,
        :forecast_string => "",
        :date_inds       => OrderedDict(d => i for (i, d) in enumerate(dates)),
        :indices         => OrderedDict(v => i for (i, v) in enumerate(vars)),
    )
    global mb_full = MeansBands(metadata, means_df, bands_dict)
end

@testset "Test that you can construct MeansBands objects and do stuff with them" begin
    # If one is empty, just return the non-empty
    @test cat(mb_empty, mb_full).means == mb_full.means
    @test cat(mb_full, mb_empty).means == mb_full.means
    @test get_class(mb_full) == :obs
    @test get_class(mb_empty) == :none
    @test get_product(mb_full) == :forecast
    @test get_product(mb_empty) == :none
    @test DSGE.get_cond_type(mb_full) == :none
    @test DSGE.get_cond_type(mb_empty) == :none
    @test DSGE.get_para(mb_full) == :full
    @test DSGE.get_para(mb_empty) == :none
    @test DSGE.n_vars_means(mb_full) == 3
    @test DSGE.n_vars_means(mb_empty) == 1
    @test DSGE.get_vars_means(mb_empty) == [:none]
    @test DSGE.get_vars_means(mb_full) == [:obs_gdp, :obs_cpi, :obs_nominalrate]
    @test DSGE.n_periods_means(mb_full) == 60
    @test DSGE.n_periods_means(mb_empty) == 1
    @test DSGE.startdate_means(mb_full) == quartertodate("2015-Q4")
    @test DSGE.startdate_means(mb_empty) == Dates.Date(0000, 1, 1)
    @test DSGE.enddate_means(mb_full) == quartertodate("2030-Q3")
    @test DSGE.enddate_means(mb_empty) == Dates.Date(0000, 1, 1)

    @test DSGE.n_vars_bands(mb_full) == 3
    @test DSGE.n_vars_bands(mb_empty) == 1
    @test DSGE.n_periods_bands(mb_full) == 60
    @test DSGE.n_periods_bands(mb_empty) == 1
    @test DSGE.startdate_bands(mb_full) == quartertodate("2015-Q4")
    @test DSGE.startdate_bands(mb_empty) == Dates.Date(0000, 1, 1)
    @test DSGE.enddate_bands(mb_full) == quartertodate("2030-Q3")
    @test DSGE.enddate_bands(mb_empty) == Dates.Date(0000, 1, 1)

    # Can only call on shockdec or irf MeansBands objects
    @test_throws AssertionError DSGE.get_shocks(mb_empty)
    @test_throws AssertionError DSGE.get_variables(mb_empty)
    @test_throws AssertionError DSGE.get_scenario_key(mb_empty)

    @test which_density_bands(mb_full) == ["90.0% LB", "80.0% LB", "70.0% LB", "60.0% LB", "50.0% LB", "50.0% UB", "60.0% UB", "70.0% UB", "80.0% UB", "90.0% UB"]
    @test which_density_bands(mb_full, uniquify = true) == ["50.0%", "60.0%", "70.0%", "80.0%", "90.0%"]
    @test which_density_bands(mb_empty) == String[]

    DSGE.prepare_meansbands_table_timeseries(mb_full, :obs_gdp)

    DSGE.parse_transform(Symbol("DSGE.loggrowthtopct_annualized"))([1.0]) == loggrowthtopct_annualized([1.0])

    # Ideally would also test prepare_means_table_shockdec...
end

@testset "Test meansbands_to_matrix works" begin
    meansbands_to_matrix(mb_full)
end

################
# Benchmarking #
################
run_benchmarks = false

if run_benchmarks
    b_cat       = @benchmark cat($mb_full, $mb_empty)
    b_to_matrix = @benchmark meansbands_to_matrix($mb_full)
    b_table     = @benchmark DSGE.prepare_meansbands_table_timeseries($mb_full, :obs_gdp)

    println("\n===== analysis/meansbands benchmark results =====")
    println("cat                              time:   ", BenchmarkTools.prettytime(median(b_cat).time),
            "   memory: ", BenchmarkTools.prettymemory(median(b_cat).memory))
    println("meansbands_to_matrix             time:   ", BenchmarkTools.prettytime(median(b_to_matrix).time),
            "   memory: ", BenchmarkTools.prettymemory(median(b_to_matrix).memory))
    println("prepare_meansbands_table_timeseries  time:   ", BenchmarkTools.prettytime(median(b_table).time),
            "   memory: ", BenchmarkTools.prettymemory(median(b_table).memory))
end
