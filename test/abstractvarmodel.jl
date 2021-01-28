# Check overloaded access functions
m = DSGEVAR(AnSchorfheide(), [:z_sh])

@testset "Access functions" begin
    # Number of anticipated policy shocks
    @test n_mon_anticipated_shocks(m) == get_setting(m, :n_mon_anticipated_shocks)
    @test n_mon_anticipated_shocks_padding(m) == get_setting(m, :n_mon_anticipated_shocks_padding)

    # Dates, indices, number of periods for each regime
    @test date_presample_start(m) == get_setting(m, :date_presample_start)
    @test date_mainsample_start(m) == get_setting(m, :date_mainsample_start)
    @test date_zlb_start(m) == get_setting(m, :date_zlb_start)
    # date_zlb_end(m) == get_setting(m, :date_zlb_end)

    @test date_presample_end(m) == Dates.lastdayofquarter(get_setting(m, :date_mainsample_start) - Dates.Month(3))
    @test date_prezlb_end(m) == Dates.lastdayofquarter(get_setting(m, :date_zlb_start) - Dates.Month(3))
    @test date_mainsample_end(m) == Dates.lastdayofquarter(get_setting(m, :date_forecast_start) - Dates.Month(3))
    @test date_conditional_end(m) == get_setting(m, :date_conditional_end)

    @test index_presample_start(m) == 1
    @test index_mainsample_start(m) == subtract_quarters(date_mainsample_start(m), date_presample_start(m)) + 1
    @test index_zlb_start(m) == subtract_quarters(date_zlb_start(m), date_presample_start(m)) + 1
    @test index_forecast_start(m) == subtract_quarters(date_forecast_start(m), date_presample_start(m)) + 1

    # Interface for fields and size
    @test DSGE.n_lags(m)        == m.lags
    @test DSGE.get_lags(m)      == m.lags
    @test DSGE.get_shocks(m)    == m.shocks
    @test DSGE.get_λ(m)         == m.λ
    @test DSGE.get_dsge(m)      == m.dsge
    @test DSGE.n_observables(m) == length(m.observables)
    @test DSGE.n_shocks(m)      == length(m.shocks)

    # Interface for accessing parameters
    @test DSGE.get_parameters(m) == m.dsge.parameters
    @test n_parameters(m)   == n_parameters(m.dsge)
    @test n_parameters_steady_state(m) == n_parameters_steady_state(m.dsge)
    @test n_parameters_free(m) == n_parameters_free(m.dsge)

    # Interface for accessing rng
    @test DSGE.get_rng(m) == m.dsge.rng

    # Interface for accessing settings dictionary
    @test DSGE.get_settings(m) == m.dsge.settings

    # Interface for accessing observables dictionary
    @test DSGE.get_observables(m) == m.observables

    # Interface for data
    @test cond_vintage(m)    == get_setting(m, :cond_vintage)
    @test cond_id(m)         == get_setting(m, :cond_id)
    @test cond_full_names(m) == get_setting(m, :cond_full_names)
    @test cond_semi_names(m) == get_setting(m, :cond_semi_names)
    @test use_population_forecast(m) == get_setting(m, :use_population_forecast)
    @test DSGE.hpfilter_population(m)     == get_setting(m, :hpfilter_population)

    # Interface for general computation settings
    @test use_parallel_workers(m)    == get_setting(m, :use_parallel_workers)

    # Interface for estimation settings
    @test reoptimize(m)          == get_setting(m, :reoptimize)
    @test calculate_hessian(m) == get_setting(m, :calculate_hessian)
    @test_throws KeyError DSGE.hessian_path(m)      == get_setting(m, :hessian_path)
    m <= Setting(:hessian_path, "")
    @test DSGE.hessian_path(m)      == get_setting(m, :hessian_path)

    @test n_hessian_test_params(m) == get_setting(m, :n_hessian_test_params)


    # Interface for Metropolis-Hastings settings
    @test n_mh_blocks(m)       == get_setting(m, :n_mh_blocks)
    @test DSGE.n_mh_param_blocks(m) == get_setting(m, :n_mh_param_blocks)
    @test n_mh_simulations(m)  == get_setting(m, :n_mh_simulations)
    @test n_mh_burn(m)         == get_setting(m, :n_mh_burn)
    @test mh_thin(m)           == get_setting(m, :mh_thin)

    # Interface for forecast settings
    @test date_forecast_start(m)   == get_setting(m, :date_forecast_start)
    @test forecast_block_size(m)   == get_setting(m, :forecast_block_size)
    @test forecast_start_block(m)  == get_setting(m, :forecast_start_block)
    @test forecast_input_file_overrides(m) == get_setting(m, :forecast_input_file_overrides)
    @test isnull(forecast_uncertainty_override(m)) # == get_setting(m, :forecast_uncertainty_override)
    @test forecast_smoother(m)     == get_setting(m, :forecast_smoother)
    @test forecast_tdist_df_val(m) == get_setting(m, :forecast_tdist_df_val)
    @test forecast_tdist_shocks(m) == get_setting(m, :forecast_tdist_shocks)
    @test forecast_zlb_value(m)    == get_setting(m, :forecast_zlb_value)
    @test impulse_response_horizons(m) == get_setting(m, :impulse_response_horizons)
    DSGE.get_forecast_input_file(m, :mode) == DSGE.get_forecast_input_file(m.dsge, :mode)

    # Interface for alternative policy settings
    m.dsge <= Setting(:alternative_policy, DSGE.default_policy())
    @test alternative_policy(m) == get_setting(m, :alternative_policy)
    delete!(m.dsge.settings, :alternative_policy)

    @test date_forecast_end(m) == Dates.lastdayofquarter(date_forecast_start(m) + Dates.Month(3 * (forecast_horizons(m)-1)))

    @test DSGE.forecast_horizons(m) == get_setting(m, :forecast_horizons)
    @test DSGE.forecast_horizons(m; cond_type = :full) == get_setting(m, :forecast_horizons) - DSGE.n_conditional_periods(m.dsge)
end

@testset "Updating parameters and settings, transforming parameters" begin
    dsge = Model1002("ss10", testing = true)
    m = DSGEVAR(dsge, [:ztil_sh])
    oldvals = map(x -> x.value, DSGE.get_parameters(m))
    newvals = oldvals + randn(n_parameters(m)) .* 1e-3
    DSGE.update!(dsge, oldvals)
    DSGE.transform_to_model_space!(m, ModelConstructors.transform_to_real_line(DSGE.get_parameters(m),
                                                                               newvals))
    DSGE.transform_to_model_space!(dsge, ModelConstructors.transform_to_real_line(DSGE.get_parameters(dsge),
                                                                                  newvals))
    @test map(x -> x.value, DSGE.get_parameters(m)) ≈ map(x -> x.value, DSGE.get_parameters(dsge))
    @test map(x -> x.value, DSGE.get_dsge(m).steady_state) ≈ map(x -> x.value, dsge.steady_state)

    DSGE.update!(dsge, oldvals)
    DSGE.update!(m, oldvals)
    @test map(x -> x.value, DSGE.get_parameters(m)) ≈ map(x -> x.value, DSGE.get_parameters(dsge))
    @test map(x -> x.value, DSGE.get_dsge(m).steady_state) ≈ map(x -> x.value, dsge.steady_state)

    newvals = oldvals + randn(n_parameters(m)) .* 1e-3
    oldpvec = dsge.parameters
    DSGE.update!(dsge, newvals)
    DSGE.update!(m, newvals)
    DSGE.update!(dsge, oldpvec)
    DSGE.update!(m, oldpvec)
    @test map(x -> x.value, DSGE.get_parameters(m)) ≈ map(x -> x.value, DSGE.get_parameters(dsge))
    @test map(x -> x.value, DSGE.get_dsge(m).steady_state) ≈ map(x -> x.value, dsge.steady_state)

    @test DSGE.get_parameters(m)[DSGE.get_dsge(m).keys[:α]].value != 0.5
    m <= parameter(:α, 0.5)
    @test DSGE.get_parameters(m)[DSGE.get_dsge(m).keys[:α]].value == 0.5

    @test DSGE.get_dsge(m).steady_state[DSGE.get_dsge(m).keys[:z_star] - n_parameters(m)].value != 0.5
    m <= SteadyStateParameter(:z_star, 0.5)
    @test DSGE.get_dsge(m).steady_state[DSGE.get_dsge(m).keys[:z_star] - n_parameters(m)].value == 0.5

    @test !haskey(DSGE.get_settings(m), :test_setting) && !haskey(DSGE.get_dsge(m).test_settings, :test_setting)
    m <= Setting(:test_setting, true)
    @test get_setting(m, :test_setting) && haskey(DSGE.get_dsge(m).test_settings, :test_setting)
    dsge = Model1002("ss10")
    m = DSGEVAR(dsge, [:ztil_sh])
    @test !haskey(DSGE.get_settings(m), :test_setting) && !haskey(DSGE.get_dsge(m).test_settings, :test_setting)
    m <= Setting(:test_setting, true)
    @test get_setting(m, :test_setting) && haskey(DSGE.get_settings(m), :test_setting)

    @test prior(m) == prior(DSGE.get_dsge(m))
end

@testset "File paths" begin
    m = DSGEVAR(AnSchorfheide("ss0", testing = true), [:z_sh])
    for folder in ["raw", "work", "figures", "tables", "log"]
        @test eval(Symbol(folder, "path"))(m, "estimate") ==
            joinpath(saveroot(m), "output_data", "an_schorfheide", "ss0", "dsgevar_an_schorfheide", "ss0", "estimate", folder)
    end
end
