test_date = Date("2000-03-31", "yyyy-mm-dd")
test_date_plus2qtr = Date("2000-09-30", "yyyy-mm-dd")
df_date = DataFrame()
df_date[!,:date] = ["2000-03-31", "2000-04-30"]
@testset "Test quarter and date utility functions" begin
    @test DSGE.prev_quarter(test_date) == Date("1999-12-31", "yyyy-mm-dd")
    @test DSGE.next_quarter(test_date) == Date("2000-06-30", "yyyy-mm-dd")
    qtr_end_ans = Vector{Date}(undef,3)
    qtr_end_ans[1] = test_date
    qtr_end_ans[2] = Date("2000-06-30", "yyyy-mm-dd")
    qtr_end_ans[3] = test_date_plus2qtr
    @test DSGE.get_quarter_ends(test_date, test_date_plus2qtr) == qtr_end_ans
    @test DSGE.quartertodate("2000-Q1") == test_date
    @test DSGE.datetoquarter(test_date) == 1
    @test DSGE.vinttodate("010331") == Date("2001-03-31", "yyyy-mm-dd")
    @test DSGE.subtract_quarters(test_date_plus2qtr, test_date) == 2
    @test DSGE.format_dates!(:date, df_date) == qtr_end_ans[1:2]
    @test DSGE.iterate_quarters(test_date, 2) == qtr_end_ans[3]
end

custom_settings = Dict{Symbol, Setting}(
    :data_vintage             => Setting(:data_vintage, "160812"),
    :cond_vintage             => Setting(:cond_vintage, "160812"),
    :cond_id                  => Setting(:cond_id, 0),
    :use_population_forecast  => Setting(:use_population_forecast, true),
    :date_forecast_start      => Setting(:date_forecast_start, DSGE.quartertodate("2016-Q3")),
    :date_conditional_end     => Setting(:date_forecast_start, DSGE.quartertodate("2016-Q3")),
    :n_anticipated_shocks     => Setting(:n_anticipated_shocks, 6))

m = Model990(custom_settings = custom_settings, testing = true)
m <= Setting(:rate_expectations_source, :ois)
@testset "Test NaN and missing data checks" begin
    @test sum(isnan.(DSGE.missing2nan(missings(3)))) == 3
    df_nan2miss = DataFrame(:a => vcat(1., NaN))
    DSGE.nan2missing!(df_nan2miss)
    @test sum(ismissing.(df_nan2miss[!,:a])) == 1

    @info "The following warnings are expected."
    df_full = load_data(m; cond_type = :full, check_empty_columns = false,
                        verbose = :none, summary_statistics = :none)
    df_semi = load_data(m; cond_type = :semi, check_empty_columns = false,
                        verbose = :none, summary_statistics = :none)

    @test_throws ErrorException load_data(m; cond_type = :full, summary_statistics = :none,
                                          verbose = :none)
    @test_throws ErrorException load_data(m; cond_type = :semi, summary_statistics = :none,
                                          verbose = :none)

    df1, df2 = DSGE.reconcile_column_names(DataFrame(a = [1, 2]), DataFrame(b = [2, 3]))
    @test sum(ismissing.(df1[!,:b])) == 2 && sum(ismissing.(df2[!,:a])) == 2

    m <= Setting(:dataroot, "")
    @test get_data_filename(m, :semi) == "data/data_cdid=00_cdvt=160812_dsid=02_vint=160812.csv"
    @test get_data_filename(m, :full) == "data/data_cdid=00_cdvt=160812_dsid=02_vint=160812.csv"
    @test get_data_filename(m, :none) == "data/data_dsid=02_vint=160812.csv"
end

nothing
