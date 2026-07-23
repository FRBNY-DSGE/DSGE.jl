path = dirname(@__FILE__)
using DSGE, BenchmarkTools, JLD2, Dates, DataFrames, OrderedCollections

m = AnSchorfheide()

# Build a minimal MeansBands in-memory (avoids stale JLD2 reference files)
function _make_io_mb(; n = 5, vars = [:obs_gdp], product = :forecast)
    dates = [Dates.lastdayofquarter(Date(2020, 1, 1) + Dates.Month(3*(i-1))) for i in 1:n]
    means_df   = DataFrame(:date => dates, (v => zeros(n) for v in vars)...)
    bands_dict = Dict{Symbol,DataFrame}(
        v => DataFrame(:date => dates, Symbol("16.0%") => zeros(n), Symbol("84.0%") => zeros(n))
        for v in vars
    )
    metadata = Dict{Symbol,Any}(
        :product         => product,
        :class           => :obs,
        :input_type      => :full,
        :cond_type       => :none,
        :para            => :full,
        :forecast_string => "",
        :date_inds       => OrderedDict(d => i for (i, d) in enumerate(dates)),
        :indices         => OrderedDict(v => i for (i, v) in enumerate(vars)),
    )
    MeansBands(metadata, means_df, bands_dict)
end

function _save_mb_tmp(mb::MeansBands)
    fn = tempname() * ".jld2"
    jldopen(fn, "w") do f; f["mb"] = mb; end
    fn
end

@testset "Testing getting MeansBands Input and Output Files" begin
    @test get_meansbands_input_file(m, :mode, :none, :histobs) ==
        joinpath(rawpath(m, "forecast"), "histobs_cond=none_para=mode_vint=" * data_vintage(m) * ".jld2")
    @test get_meansbands_input_file("a", ["b"], :mode, :none, :histobs, fileformat = :jld2) ==
        joinpath("a", "histobs_b_cond=none_para=mode.jld2")

    @test get_meansbands_output_file(m, :mode, :none, :histobs) ==
        joinpath(workpath(m, "forecast"), "mbhistobs_cond=none_para=mode_vint=" * data_vintage(m) * ".jld2")
    @test get_meansbands_output_file("a", ["b"], :mode, :none, :histobs, fileformat = "jld2") ==
        joinpath("a", "mbhistobs_b_cond=none_para=mode.jld2")

    run(`rm -r a`)

    # Test read_mb(fn): verify round-trip through JLD2
    mb_orig = _make_io_mb()
    fn = _save_mb_tmp(mb_orig)
    mb_loaded = read_mb(fn)
    @test mb_loaded isa MeansBands
    @test mb_loaded.means == mb_orig.means

    # Test read_mb(fn1, fn2) with fn1 == fn2: uses metadata from fn1, means from fn2
    mb_loaded2 = read_mb(fn, fn)
    @test mb_loaded2 isa MeansBands
    @test mb_loaded2.means == mb_orig.means

    # Test read_bdd_and_unbdd_mb: empty filenames throw
    @test_throws AssertionError read_bdd_and_unbdd_mb("", "")

    # Test read_bdd_and_unbdd_mb: result has unbdd means + bdd bands
    mb_bdd   = _make_io_mb(product = :bddforecast)
    mb_unbdd = _make_io_mb(product = :forecast)
    fn_bdd   = _save_mb_tmp(mb_bdd)
    fn_unbdd = _save_mb_tmp(mb_unbdd)
    mb_merged = read_bdd_and_unbdd_mb(fn_bdd, fn_unbdd)
    @test mb_merged isa MeansBands
    @test mb_merged.means == mb_unbdd.means
    @test mb_merged.bands == mb_bdd.bands
end

@testset "Testing add_requisite_output_vars_meansbands" begin
    @test DSGE.add_requisite_output_vars_meansbands([:histobs, :shockdecpseudo]) ==
        [:histobs, :shockdecpseudo, :dettrendpseudo, :trendpseudo, :histforecastpseudo]
    @test DSGE.add_requisite_output_vars_meansbands([:histobs, :shockdecobs]) ==
        [:histobs, :shockdecobs, :dettrendobs, :trendobs, :histforecastobs]
end

################
# Benchmarking #
################
run_benchmarks = false

if run_benchmarks
    mb_bench     = _make_io_mb(n = 60, vars = [:obs_gdp, :obs_cpi, :obs_nominalrate])
    mb_bdd_bench = _make_io_mb(n = 60, vars = [:obs_gdp, :obs_cpi, :obs_nominalrate], product = :bddforecast)
    mb_ubd_bench = _make_io_mb(n = 60, vars = [:obs_gdp, :obs_cpi, :obs_nominalrate], product = :forecast)
    fn_bench     = _save_mb_tmp(mb_bench)
    fn_bdd_bench = _save_mb_tmp(mb_bdd_bench)
    fn_ubd_bench = _save_mb_tmp(mb_ubd_bench)

    b_read_mb = @benchmark read_mb($fn_bench)
    b_bdd     = @benchmark read_bdd_and_unbdd_mb($fn_bdd_bench, $fn_ubd_bench)

    println("\n===== analysis/io benchmark results =====")
    println("read_mb                time:   ", BenchmarkTools.prettytime(median(b_read_mb).time),
            "   memory: ", BenchmarkTools.prettymemory(median(b_read_mb).memory))
    println("read_bdd_and_unbdd_mb  time:   ", BenchmarkTools.prettytime(median(b_bdd).time),
            "   memory: ", BenchmarkTools.prettymemory(median(b_bdd).memory))
end
