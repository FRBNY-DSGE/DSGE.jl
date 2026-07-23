using DSGE, Dates, DataFrames, BenchmarkTools
path = dirname(@__FILE__)

# Construct a synthetic MeansBands in-memory to avoid stale JLD2 reference files.
# 60 quarterly periods starting 2015-Q4; two observables.
let n = 60, vars = [:obs_gdp, :obs_cpi]
    dates = [Dates.lastdayofquarter(Date(2015, 10, 1) + Dates.Month(3*(i-1))) for i in 1:n]
    means_df = DataFrame(:date => dates, (v => zeros(n) for v in vars)...)
    bands_dict = Dict{Symbol,DataFrame}(
        v => DataFrame(:date => dates, Symbol("16.0%") => zeros(n), Symbol("84.0%") => zeros(n))
        for v in vars
    )
    metadata = Dict{Symbol,Any}(
        :product    => :forecast,
        :class      => :obs,
        :input_type => :full,
        :cond_type  => :none,
        :date_inds  => Dict(d => i for (i, d) in enumerate(dates)),
        :indices    => Dict(v => i for (i, v) in enumerate(vars)),
    )
    global mb_full = MeansBands(metadata, means_df, bands_dict)
end

# create_q4q4_mb requires a 4q product — make a copy with the right product
mb_4q = deepcopy(mb_full)
mb_4q.metadata[:product] = :forecast4q

mb_q4q4 = create_q4q4_mb(mb_4q)

@testset "create_q4q4_mb" begin
    # Error on non-4q product
    @test_throws ErrorException create_q4q4_mb(mb_full)

    # Product name updated correctly (:forecast4q -> :forecastq4q4)
    @test mb_q4q4.metadata[:product] == :forecastq4q4

    # All dates in metadata are Q4
    @test all(Dates.quarterofyear(d) == 4 for d in keys(mb_q4q4.metadata[:date_inds]))

    # All rows in means are Q4
    @test all(Dates.quarterofyear(d) == 4 for d in mb_q4q4.means[!, :date])

    # Fewer rows than original (Q4 only is a subset)
    @test size(mb_q4q4.means, 1) < size(mb_4q.means, 1)

    # Bands also filtered to Q4
    for var in keys(mb_q4q4.bands)
        @test all(Dates.quarterofyear(d) == 4 for d in mb_q4q4.bands[var][!, :date])
    end

    # All valid product types pass through without error
    for prod in [:hist4q, :histforecast4q, :bddforecast4q, :bddhistforecast4q]
        mb_test = deepcopy(mb_full)
        mb_test.metadata[:product] = prod
        result = create_q4q4_mb(mb_test)
        @test all(Dates.quarterofyear(d) == 4 for d in result.means[!, :date])
    end
end

################
# Benchmarking #
################
run_benchmarks = false

if run_benchmarks
    b_q4q4 = @benchmark create_q4q4_mb($mb_4q)

    println("\n===== create_q4q4_mb benchmark results =====")
    println(rpad("create_q4q4_mb", 18), " time: ", rpad(BenchmarkTools.prettytime(median(b_q4q4).time), 12),
            "memory: ", BenchmarkTools.prettymemory(median(b_q4q4).memory))
end
