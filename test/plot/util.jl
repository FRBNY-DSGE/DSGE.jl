quarter_date_to_number_out = DSGE.quarter_date_to_number(DSGE.quartertodate("2019-Q1"))
quarter_date_to_number_ans = 2019.00
quarter_number_to_date_out = DSGE.quarter_number_to_date(quarter_date_to_number_ans)
quarter_number_to_date_ans = DSGE.quartertodate("2019-Q1")
dates = [DSGE.quartertodate("2019-Q1"),DSGE.quartertodate("2019-Q2"),DSGE.quartertodate("2019-Q3"),DSGE.quartertodate("2019-Q4")]
date_to_float_out = date_to_float(dates, 2)
date_to_float_ans = [2019.25 2019.5 2019.75 2020.;
          2019.25 2019.5 2019.75 2020.]
m = AnSchorfheide(testing = true)
transform_list = [DSGE.loggrowthtopct_annualized_percapita, DSGE.loggrowthtopct_annualized,
                  DSGE.loggrowthtopct, DSGE.loggrowthtopct_percapita,
                  DSGE.logleveltopct_annualized_percapita, DSGE.logleveltopct_annualized,
                  DSGE.quartertoannual, DSGE.loggrowthtopct, identity]
untrans_ans = ["Q/Q Log Growth Rate", "Q/Q Log Growth Rate", "Q/Q Log Growth Rate",
               "Q/Q Log Growth Rate", "Log Level",  "Log Level", "Percent Q/Q",
               "Q/Q Log Growth Rate", ""]
q4_ans = ["Percent 4Q Growth", "Percent 4Q Growth", "Percent 4Q Growth",
          "Percent 4Q Growth", "Percent 4Q Growth",  "Percent 4Q Growth",
          "Percent Annualized", "Percent 4Q Growth", ""]
other_ans = ["Percent Q/Q Annualized", "Percent Q/Q Annualized", "Percent Q/Q Annualized",
          "Percent Q/Q Annualized", "Percent Q/Q Annualized",  "Percent Q/Q Annualized",
          "Percent Annualized", "Percent Q/Q Annualized", ""]

@testset "Check utility functions for plotting" begin
    @test "Real GDP Growth" == DSGE.describe_series(m, :obs_gdp, :obs)
    @test "Output Growth" == DSGE.describe_series(m, :y_t, :pseudo)
    @test "π_t" == DSGE.describe_series(m, :π_t, :states)
    @test "pi_t" == DSGE.describe_series(m, :π_t, :states, detexify = true)
    @test "z" == DSGE.describe_series(m, :z_sh, :shocks)
    for (i,tf) in enumerate(transform_list)
        m.observable_mappings[:obs_gdp].rev_transform = tf
        m.pseudo_observable_mappings[:y_t].rev_transform = tf
        @test untrans_ans[i] == DSGE.series_ylabel(m, :obs_gdp, :obs, untrans = true)
        @test untrans_ans[i] == DSGE.series_ylabel(m, :y_t, :pseudo, untrans = true)
        @test q4_ans[i] == DSGE.series_ylabel(m, :obs_gdp, :obs, fourquarter = true)
        @test q4_ans[i] == DSGE.series_ylabel(m, :y_t, :pseudo, fourquarter = true)
        @test other_ans[i] == DSGE.series_ylabel(m, :obs_gdp, :obs)
        @test other_ans[i] == DSGE.series_ylabel(m, :y_t, :pseudo)
    end
    m.observable_mappings[:obs_gdp].rev_transform = sum
    @test_throws ErrorException DSGE.series_ylabel(m, :obs_gdp, :obs)
    @test_throws ErrorException DSGE.series_ylabel(m, :obs_gdp, :test)
    @test "Standard Deviations" == DSGE.series_ylabel(m, :test, :stdshocks)
    @test isempty(DSGE.series_ylabel(m, :test, :states))
    @test isempty(DSGE.series_ylabel(m, :test, :shocks))
    @test quarter_date_to_number_out == quarter_date_to_number_ans
    @test quarter_number_to_date_out == quarter_number_to_date_ans
    @test date_to_float_out == date_to_float_ans
end

nothing
