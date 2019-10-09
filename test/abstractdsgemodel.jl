using DSGE, ModelConstructors, Dates, Optim, Test, Nullables, HDF5
path = dirname(@__FILE__)
m = AnSchorfheide()
m <= Setting(:date_zlb_end, get_setting(m, :date_zlb_start))
m <= Setting(:hessian_path, "")
# Dates, indices, number of periods for each regime
@testset "Test setting access functions for AbstractDSGEModel objects" begin
    n_anticipated_shocks(m) == get_setting(m, :n_anticipated_shocks)
    n_anticipated_shocks_padding(m) == get_setting(m, :n_anticipated_shocks_padding)

    @test DSGE.date_presample_start(m) == get_setting(m, :date_presample_start)
    @test DSGE.date_mainsample_start(m) == get_setting(m, :date_mainsample_start)
    @test DSGE.date_zlb_start(m) == get_setting(m, :date_zlb_start)
    @test DSGE.date_zlb_end(m) == get_setting(m, :date_zlb_end)

    @test DSGE.date_presample_end(m) == Dates.lastdayofquarter(get_setting(m, :date_mainsample_start) - Dates.Month(3))
    @test DSGE.date_prezlb_end(m) == Dates.lastdayofquarter(get_setting(m, :date_zlb_start) - Dates.Month(3))
    @test DSGE.date_mainsample_end(m) == Dates.lastdayofquarter(get_setting(m, :date_forecast_start) - Dates.Month(3))
    @test DSGE.date_conditional_end(m) == get_setting(m, :date_conditional_end)

    @test DSGE.index_presample_start(m) == 1
    @test DSGE.index_mainsample_start(m) == subtract_quarters(DSGE.date_mainsample_start(m), DSGE.date_presample_start(m)) + 1
    @test DSGE.index_zlb_start(m) == subtract_quarters(DSGE.date_zlb_start(m), DSGE.date_presample_start(m)) + 1
    @test DSGE.index_forecast_start(m) == subtract_quarters(DSGE.date_forecast_start(m), DSGE.date_presample_start(m)) + 1

    @test DSGE.index_shockdec_start(m) == subtract_quarters(DSGE.date_shockdec_start(m), DSGE.date_mainsample_start(m)) + 1

    @test DSGE.index_shockdec_end(m) == subtract_quarters(DSGE.date_shockdec_end(m), DSGE.date_mainsample_start(m)) + 1

    @test n_presample_periods(m)   == subtract_quarters(DSGE.date_mainsample_start(m), DSGE.date_presample_start(m))
    @test n_prezlb_periods(m)      == subtract_quarters(DSGE.date_zlb_start(m), DSGE.date_mainsample_start(m))
    @test n_zlb_periods(m)         == subtract_quarters(DSGE.date_forecast_start(m), DSGE.date_zlb_start(m))
    @test n_mainsample_periods(m)  == subtract_quarters(DSGE.date_forecast_start(m), DSGE.date_mainsample_start(m))
    @test n_conditional_periods(m) == subtract_quarters(DSGE.date_conditional_end(m), DSGE.date_mainsample_end(m))

    @test inds_presample_periods(m)  == collect(DSGE.index_presample_start(m):(DSGE.index_mainsample_start(m)-1))
    @test inds_prezlb_periods(m)     == collect(DSGE.index_mainsample_start(m):(DSGE.index_zlb_start(m)-1))
    @test inds_zlb_periods(m)        == collect(DSGE.index_zlb_start(m):(DSGE.index_forecast_start(m)-1))
    @test inds_mainsample_periods(m) == collect(DSGE.index_mainsample_start(m):(DSGE.index_forecast_start(m)-1))

    @test cond_vintage(m)    == get_setting(m, :cond_vintage)
    @test cond_id(m)         == get_setting(m, :cond_id)
    @test cond_full_names(m) == get_setting(m, :cond_full_names)
    @test cond_semi_names(m) == get_setting(m, :cond_semi_names)
    @test use_population_forecast(m) == get_setting(m, :use_population_forecast)
    @test DSGE.hpfilter_population(m)     == get_setting(m, :hpfilter_population)

    # Interface for general computation settings
    @test DSGE.use_parallel_workers(m)    == get_setting(m, :use_parallel_workers)

    # Interface for estimation settings
    @test DSGE.reoptimize(m)          == get_setting(m, :reoptimize)
    @test DSGE.calculate_hessian(m) == get_setting(m, :calculate_hessian)
    @test DSGE.hessian_path(m)      == get_setting(m, :hessian_path)
    @test DSGE.n_hessian_test_params(m) == get_setting(m, :n_hessian_test_params)

    # Interface for Metropolis-Hastings settings
    @test DSGE.n_mh_blocks(m)       ==  get_setting(m, :n_mh_blocks)
    @test DSGE.n_mh_param_blocks(m) ==  get_setting(m, :n_mh_param_blocks)
    @test DSGE.n_mh_simulations(m)  ==  get_setting(m, :n_mh_simulations)
    @test DSGE.n_mh_burn(m)         ==  get_setting(m, :n_mh_burn)
    @test DSGE.mh_thin(m)           ==  get_setting(m, :mh_thin)

    # Interface for forecast settings
    @test DSGE.date_forecast_start(m)   == get_setting(m, :date_forecast_start)
    @test DSGE.forecast_block_size(m)   == get_setting(m, :forecast_block_size)
    @test DSGE.forecast_start_block(m)  == get_setting(m, :forecast_start_block)
    @test DSGE.forecast_input_file_overrides(m) == get_setting(m, :forecast_input_file_overrides)
    @test isnull(DSGE.forecast_uncertainty_override(m)) && isnull(get_setting(m, :forecast_uncertainty_override))
    @test DSGE.forecast_smoother(m)     == get_setting(m, :forecast_smoother)
    @test DSGE.forecast_tdist_df_val(m) == get_setting(m, :forecast_tdist_df_val)
    @test DSGE.forecast_tdist_shocks(m) == get_setting(m, :forecast_tdist_shocks)
    @test DSGE.forecast_zlb_value(m)    == get_setting(m, :forecast_zlb_value)
    @test DSGE.impulse_response_horizons(m) == get_setting(m, :impulse_response_horizons)
    @test DSGE.n_shockdec_periods(m)    == DSGE.index_shockdec_end(m) - DSGE.index_shockdec_start(m) + 1

    # Interface for alternative policy settings
    @test DSGE.alternative_policy(m) == get_setting(m, :alternative_policy)
end

@testset "Test other auxiliary setting functions for AbstractDSGEModel objects" begin
    m <= Setting(:population_mnemonic, :test)
    @test get(DSGE.parse_population_mnemonic(m)[1]) == :test
    m <= Setting(:population_mnemonic, Nullable())
    println(DSGE.parse_population_mnemonic(m))
    @test DSGE.parse_population_mnemonic(m)[1].hasvalue == false
    m1 = Model1002("ss10")
    @test @test_matrix_approx_eq inds_states_no_ant(m) collect(1:8)
    @test @test_matrix_approx_eq inds_states_no_ant(m1) vcat(collect(1:62), collect(69:84))
    @test @test_matrix_approx_eq inds_shocks_no_ant(m) collect(1:3)
    @test @test_matrix_approx_eq inds_shocks_no_ant(m1) collect(1:18)
    @test @test_matrix_approx_eq inds_obs_no_ant(m) collect(1:3)
    @test @test_matrix_approx_eq inds_obs_no_ant(m1) collect(1:13)


    m <= Setting(:date_forecast_start, DSGE.quartertodate("2019-Q3"))
    m <= Setting(:forecast_horizons, 10)
    exp_date_forecast_end = DSGE.lastdayofquarter(DSGE.quartertodate("2019-Q3") + Dates.Month(3 * 9))
    @test forecast_horizons(m) == 10
    m <= Setting(:date_conditional_end, DSGE.date_mainsample_end(m) + Dates.Month(3))
    @test forecast_horizons(m; cond_type = :full) == 9
    @test forecast_horizons(m; cond_type = :semi) == 9
    @test DSGE.date_forecast_end(m) == exp_date_forecast_end
    @test DSGE.date_shockdec_start(m) == date_mainsample_start(m)
    m <= Setting(:shockdec_startdate, nothing)
    @test DSGE.date_shockdec_start(m) == DSGE.date_mainsample_start(m)
    m <= Setting(:shockdec_startdate, Nullable{Date}(DSGE.quartertodate("2019-Q3")))
    @test DSGE.date_shockdec_start(m) == DSGE.quartertodate("2019-Q3")
    m <= Setting(:shockdec_enddate, Nullable{Date}())
    @test DSGE.date_shockdec_end(m) == DSGE.date_forecast_end(m)
    m <= Setting(:shockdec_enddate, Nullable{Date}(DSGE.quartertodate("2019-Q3")))
    @test DSGE.date_shockdec_end(m) == DSGE.quartertodate("2019-Q3")

    @test @test_matrix_approx_eq load_parameters_from_file(m, "$path/reference/load_params.h5") h5read("$path/reference/load_params.h5", "params")
    @test_throws ErrorException load_parameters_from_file(m, "$path/reference/smc.h5")
    @test_throws ErrorException load_parameters_from_file(m, "$path/reference/smc.jld2")
    @test_throws AssertionError load_parameters_from_file(m1, "$path/reference/load_params.h5")
    @test_throws AssertionError load_parameters_from_file(m, "$path/reference/load_params_float32.h5")

    m.parameters[1].value = 1.
    specify_mode!(m, "reference/load_params.h5")
    out_params = map(x -> x.value, m.parameters)
    @test @test_matrix_approx_eq out_params load_parameters_from_file(m, "$path/reference/load_params.h5")

    specify_hessian!(m, "reference/hessian.h5")
    @test get_setting(m, :calculate_hessian) == false
    @test get_setting(m, :hessian_path) == "$path/reference/hessian.h5"

    para = deepcopy(map(x -> x.value, m.parameters))
    m.parameters[1].value = 1.
    DSGE.update!(m, para)
    @test m.parameters[1].value != 1.
    @test @test_matrix_approx_eq map(x->x.value,m.parameters) para
end

nothing
