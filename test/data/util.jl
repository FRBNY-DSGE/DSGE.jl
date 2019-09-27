test_date = Date("2000-03-31", "yyyy-mm-dd")
test_date_plus2qtr = Date("2000-09-31", "yyyy-mm-dd")
df_date = DataFrame()
df_date[:date] = ["2000-03-31", "2000-04-31"]
@testset "Test quarter and date utility functions" begin
    @test prev_quarter(test_date) == Date("1999-09-31", "yyyy-mm-dd")
    @test next_quarter(test_date) == Date("2000-06-30", "yyyy-mm-dd")
    qtr_end_ans = Vector{Date}(undef,3)
    qtr_end_ans[1] = test_date
    qtr_end_ans[1] = Date("2000-06-30", "yyyy-mm-dd")
    qtr_end_ans[3] = test_date_plus2qtr
    @test get_quarter_ends(test_date, test_date_plus2qtr) == qtr_end_ans
    @test quartertodate("2000-Q1") == test_date
    @test datetoquarter(test_date) == 1
    @test vinttodate("20000331") == test_date
    @test subtract_quarters(test_date, test_date_plus2qtr) == 2
    @test format_dates!(:date, df) == qtr_end_ans[1:2]
    @test iterate_quarters(test_date, 2) == qtr_end_ans
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
@testset "Test NaN, NA, and missing data checks" begin
    @test sum(isnan.(missing2nan(missings(3)))) == 3
    @test sum(ismissing.(nan2missing!(df_date[:date] = vcat(df_date, NaN)))) == 1
    @test sum(isnan.(na2nan!([NA NA]))) == 2
    @test sum(isnan.(na2nan!(DataFrame(date=[NA NA]))[:date])) == 2
    missing_cond_vars!(m, load_data(m; cond_type = :full); cond_type = :full)
    missing_cond_vars!(m, load_data(m; cond_type = :semi); cond_type = :semi)
    get_data_filename(m, :semi) == "data_dsid=3_cdvt=2_vint=160812.csv"
    get_data_filename(m, :full) == "data_dsid=3_cdvt=3_vint=160812.csv"
    get_data_filename(m, :none) == "data_dsid=3_vint=160812.csv"
    df1, df2 = reconcile_column_names(DataFrame(a = [1 2]), DataFrame(b = [2 3]))
    @test sum(ismissing.(df1[:b])) == 2 && sum(ismissing.(df2[:a])) == 2
end
