using DSGE, Test, ModelConstructors, BenchmarkTools

function make_model(year::Int)
    m = AnSchorfheide()
    m <= Setting(:saveroot, tempdir())
    m <= Setting(:data_vintage, "160101")
    m <= Setting(:date_forecast_start, DSGE.quartertodate("$year-Q1"))
    m <= Setting(:date_conditional_end, DSGE.quartertodate("$year-Q1"))
    return m
end

m_new  = make_model(2016)
m_old  = make_model(2014)   # different forecast_start → the :release component
m_same = make_model(2016)   # same forecast_start as m_new

@testset "get_decomp_filename has the new__old structure" begin
    fn = DSGE.get_decomp_filename(m_new, m_old, :full, :none, :none, :decomptotal, :obs)
    @test endswith(fn, ".jld2")
    @test occursin("__", basename(fn))                 # new__old separator
    @test occursin(spec(m_old), basename(fn))          # old model's spec is embedded
    @test occursin("decomptotalobs", basename(fn))     # product*class from the new side

    # Distinct (product, class) → distinct filenames
    fn2 = DSGE.get_decomp_filename(m_new, m_old, :full, :none, :none, :decomptrend, :pseudo)
    @test fn != fn2

    # forecast_string_new changes the new-side filename
    fn3 = DSGE.get_decomp_filename(m_new, m_old, :full, :none, :none, :decomptotal, :obs;
                                   forecast_string_new = "alt")
    @test fn != fn3
end

@testset "get_decomp_output_files spans comps × classes (with :release)" begin
    files = DSGE.get_decomp_output_files(m_new, m_old, :full, :none, :none, [:obs, :pseudo])
    # forecast starts differ → 9 comps (includes :release) × 2 classes
    @test length(files) == 9 * 2
    for v in (:decompshockdecobs, :decomptrendpseudo, :decomptotalobs, :decompreleaseobs)
        @test haskey(files, v)
    end
    @test all(endswith(f, ".jld2") for f in values(files))
end

@testset "get_decomp_output_files drops :release when forecast starts match" begin
    files = DSGE.get_decomp_output_files(m_new, m_same, :full, :none, :none, [:obs, :pseudo])
    @test length(files) == 8 * 2
    @test !haskey(files, :decompreleaseobs)
end

@testset "get_decomp_output_files: model_decomp adds the :model component" begin
    files = DSGE.get_decomp_output_files(m_new, m_old, :full, :none, :none, [:obs, :pseudo];
                                         model_decomp = true)
    @test length(files) == 10 * 2
    @test haskey(files, :decompmodelobs)
end

@testset "get_decomp_mean_file uses workpath, not rawpath" begin
    meanf = DSGE.get_decomp_mean_file(m_new, m_old, :full, :none, :none, :obs)
    @test endswith(meanf, ".jld2")
    @test occursin("work", meanf)
    rawf = DSGE.get_decomp_filename(m_new, m_old, :full, :none, :none, :decomp, :obs;
                                    pathfcn = DSGE.rawpath)
    @test meanf != rawf
end

################
# Benchmarking #
################
run_benchmarks = false

if run_benchmarks
    b_fname = @benchmark DSGE.get_decomp_filename($m_new, $m_old, :full, :none, :none, :decomptotal, :obs)
    b_files = @benchmark DSGE.get_decomp_output_files($m_new, $m_old, :full, :none, :none, [:obs, :pseudo])

    println("\n===== decomp io benchmark results =====")
    for (name, b) in [("get_decomp_filename",     b_fname),
                      ("get_decomp_output_files", b_files)]
        println(rpad(name, 26), " time: ", rpad(BenchmarkTools.prettytime(median(b).time), 12),
                "memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end
