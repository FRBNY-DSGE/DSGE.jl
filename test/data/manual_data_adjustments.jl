using DSGE, Test, ModelConstructors, DataFrames, Dates, BenchmarkTools

N = 6
m = AnSchorfheide()
m <= Setting(:n_mon_anticipated_shocks, N)
m <= Setting(:date_presample_start, Date(2019, 12, 31))

# Quarterly (end-of-quarter) dates spanning the pgap/ygap init and ZLB windows
# the function keys off of (2020-06-30 init; 2020-12-31..2021-12-31 ZLB).
dates = [Date(2020, 3, 31), Date(2020, 6, 30), Date(2020, 9, 30), Date(2020, 12, 31),
         Date(2021, 3, 31), Date(2021, 6, 30), Date(2021, 9, 30), Date(2021, 12, 31)]

function make_df()
    df = DataFrame(:date => copy(dates))
    df[!, :obs_nominalrate] = fill(1.0, length(dates))
    # The source's expected-rate masking (manual_data_adjustments.jl:53-56) maps
    # submatrix columns back with `k[2] + minimum(exp_inds) - 1`, which only works
    # if the obs_exp_nominalrate columns are CONTIGUOUS. Real data groups them, so
    # add all obs_nominalrate{i} first, then all obs_exp_nominalrate{i}.
    for i in 1:N
        df[!, "obs_nominalrate$i"] = fill(1.0, length(dates))
    end
    for i in 1:N
        df[!, "obs_exp_nominalrate$i"] = fill(0.5, length(dates))    # > 0.033
    end
    allowmissing!(df, ["obs_exp_nominalrate$i" for i in 1:N])
    return df
end

# fcast_date < 2022-04-01 avoids the (broken) mon_ant_ait_shocks branch; cond_type
# = :none and add_22Q1_ffr = false (default) keep us on the main path.
out = DSGE.post_covid_data_mods(m, make_df(), :none, DataFrame(); fcast_date = Date(2021, 1, 1))

@testset "post_covid_data_mods preserves rows when all are in-sample" begin
    @test nrow(out) == length(dates)
    @test out[!, :date] == dates
end

@testset "post_covid_data_mods adds and seeds obs_pgap / obs_ygap" begin
    @test "obs_pgap" in names(out)
    @test "obs_ygap" in names(out)
    # Seeded at the init date (2020-06-30, row 2); NaN elsewhere
    @test out[2, :obs_pgap] == -0.125
    @test out[2, :obs_ygap] == -12.0
    @test isnan(out[1, :obs_pgap])
    @test isnan(out[1, :obs_ygap])
end

@testset "post_covid_data_mods NaNs the ZLB nominal-rate window" begin
    # ZLB runs 2020-12-31..2021-12-31 (rows 4:8)
    @test all(isnan, out[4:8, :obs_nominalrate])
    @test all(out[1:3, :obs_nominalrate] .== 1.0)          # pre-ZLB untouched
    for i in 1:N
        @test all(isnan, out[4:8, "obs_nominalrate$i"])
    end
end

@testset "post_covid_data_mods NaNs the final anticipated-rate obs" begin
    # Last row's anticipated nominal-rate obs are zeroed out (cond_type = :none)
    for i in 1:N
        @test isnan(out[end, "obs_nominalrate$i"])
    end
end

@testset "post_covid_data_mods masks the expected-rate history" begin
    # keep_Q1_spd path sets rows 1:end-3 of each obs_exp_nominalrate column missing
    for i in 1:N
        @test all(ismissing, out[1:end-3, "obs_exp_nominalrate$i"])
    end
end

@testset "post_covid_data_mods drops out-of-sample (pre-presample) rows" begin
    m2 = AnSchorfheide()
    m2 <= Setting(:n_mon_anticipated_shocks, N)
    # Drop only the first row (2020-03-31); keep 2020-06-30, which the source
    # requires as the pgap/ygap init date (see the BUG testset below).
    m2 <= Setting(:date_presample_start, Date(2020, 6, 30))
    out2 = DSGE.post_covid_data_mods(m2, make_df(), :none, DataFrame(); fcast_date = Date(2021, 1, 1))
    @test nrow(out2) == length(dates) - 1
    @test minimum(out2[!, :date]) == Date(2020, 6, 30)
end

@testset "ss104_estimation returns the dataframe unchanged" begin
    df = make_df()
    @test DSGE.ss104_estimation(df) === df
end

####################
# Known source bugs #
####################
# The following are flagged @test_broken: they document real defects in
# src/data/manual_data_adjustments.jl. Each will flip to an "Unexpected Pass"
# (erroring the testset) once the source is fixed, prompting these to be updated.

@testset "BUG: mon_ant_ait_shocks branch references an undefined global" begin
    # When spd_expect_data && fcast_date >= 2022-04-01 && cond_type == :none, the
    # function hits `for i in mon_ant_ait_shocks` (manual_data_adjustments.jl:61),
    # but `mon_ant_ait_shocks` is defined nowhere in the source, so the call throws
    # an UndefVarError. Likely meant to be a kwarg or `expected_ffr`.
    @test_broken (DSGE.post_covid_data_mods(m, make_df(), :none, DataFrame();
                                            fcast_date = Date(2022, 6, 1)); true)
end

@testset "BUG: pgap/ygap init date is not nothing-guarded" begin
    # ind_init = findfirst(date == 2020-06-30) is used unguarded at line 94, unlike
    # start_ind which is guarded with `if !isnothing(...)` at line 98. If the data
    # doesn't contain 2020-06-30, ind_init is nothing and df[nothing, :obs_pgap]
    # throws. Feed a frame whose presample start excludes the init date.
    mbug = AnSchorfheide()
    mbug <= Setting(:n_mon_anticipated_shocks, N)
    mbug <= Setting(:date_presample_start, Date(2020, 9, 30))   # drops 2020-06-30
    @test_broken (DSGE.post_covid_data_mods(mbug, make_df(), :none, DataFrame();
                                            fcast_date = Date(2021, 1, 1)); true)
end

@testset "BUG: exported name post_covid_data_mods! does not exist" begin
    # DSGE.jl exports `post_covid_data_mods!` (with a bang), but the function is
    # defined as `post_covid_data_mods` (no bang), so the exported binding is
    # undefined.
    @test isdefined(DSGE, :post_covid_data_mods)                      # the real one
    @test_broken isdefined(DSGE, Symbol("post_covid_data_mods!"))     # the exported name
end

################
# Benchmarking #
################
run_benchmarks = false

if run_benchmarks
    b_mods = @benchmark DSGE.post_covid_data_mods($m, make_df(), :none, DataFrame();
                                                  fcast_date = Date(2021, 1, 1))
    b_ss104 = @benchmark DSGE.ss104_estimation(make_df())

    println("\n===== manual_data_adjustments benchmark results =====")
    for (name, b) in [("post_covid_data_mods", b_mods),
                      ("ss104_estimation",     b_ss104)]
        println(rpad(name, 24), " time: ", rpad(BenchmarkTools.prettytime(median(b).time), 12),
                "memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end
