using DSGE, Test, ModelConstructors, DataFrames, Dates, BenchmarkTools

function make_model(year::Int)
    m = AnSchorfheide()
    m <= Setting(:date_forecast_start, DSGE.quartertodate("$year-Q1"))
    m <= Setting(:date_conditional_end, DSGE.quartertodate("$year-Q1"))
    m <= Setting(:forecast_horizons, 16)
    return m
end

m_new = make_model(2016)
m_old = make_model(2014)

# Quantities the function should reproduce, derived independently from accessors.
T0     = DSGE.n_presample_periods(m_new)
T      = DSGE.n_mainsample_periods(m_new)
k      = DSGE.subtract_quarters(date_forecast_start(m_new), date_forecast_start(m_old))
T1_new = DSGE.n_conditional_periods(m_new)
T1_old = DSGE.n_conditional_periods(m_old)
H_exp  = DSGE.subtract_quarters(DSGE.date_forecast_end(m_old), DSGE.date_mainsample_end(m_new))

mkdf(n) = DataFrame(date = collect(1:n))   # only row count matters to the function

@testset "decomposition_periods: presample equal, k is the forecast-start gap" begin
    @test DSGE.n_presample_periods(m_old) == T0
    @test k == 8
end

@testset "decomposition_periods returns (T, k, H) for each cond combo" begin
    for (cond_new, cond_old) in [(:none, :none), (:none, :full), (:full, :none), (:full, :full)]
        df_new = mkdf(T0 + T + (cond_new == :full ? T1_new : 0))
        df_old = mkdf(T0 + T - k + (cond_old == :full ? T1_old : 0))
        Tr, kr, Hr = DSGE.decomposition_periods(m_new, m_old, df_new, df_old, cond_new, cond_old)
        @test Tr == T
        @test kr == k
        @test Hr == H_exp     # H does not depend on conditioning
    end
end

@testset "decomposition_periods asserts on wrong dataframe sizes" begin
    good_old = mkdf(T0 + T - k)
    # df_new too long
    @test_throws AssertionError DSGE.decomposition_periods(m_new, m_old, mkdf(T0 + T + 1), good_old, :none, :none)
    # df_old too long
    @test_throws AssertionError DSGE.decomposition_periods(m_new, m_old, mkdf(T0 + T), mkdf(T0 + T - k + 1), :none, :none)
end

@testset "decomposition_periods asserts the presample lengths match" begin
    m_bad = make_model(2014)
    m_bad <= Setting(:date_presample_start, date_presample_start(m_bad) - Dates.Month(3))
    @test_throws AssertionError DSGE.decomposition_periods(m_new, m_bad, mkdf(T0 + T), mkdf(T0 + T - k), :none, :none)
end

################
# Benchmarking #
################
run_benchmarks = false

if run_benchmarks
    df_new = mkdf(T0 + T)
    df_old = mkdf(T0 + T - k)
    b = @benchmark DSGE.decomposition_periods($m_new, $m_old, $df_new, $df_old, :none, :none)
    println("\n===== decomposition_periods benchmark =====")
    println(rpad("decomposition_periods", 24), " time: ",
            rpad(BenchmarkTools.prettytime(median(b).time), 12),
            "memory: ", BenchmarkTools.prettymemory(median(b).memory))
end
