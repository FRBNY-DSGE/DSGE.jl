m = AnSchorfheide()
DSGE.default_settings!(m)

@testset "Check default settings" begin
    # I/O File locations
    @test get_setting(m, :saveroot) == normpath(joinpath(dirname(@__FILE__), "..","save"))
    @test get_setting(m, :dataroot) == normpath(joinpath(dirname(@__FILE__), "..","save","input_data"))

    # Data settings for released and conditional data. Default behavior is to set vintage
    # of data to today's date.
    @test get_setting(m, :data_vintage) == Dates.format(now(), DSGE.DSGE_DATE_FORMAT)
    @test get_setting(m, :data_id) == 3
    @test get_setting(m, :cond_vintage) == Dates.format(now(), DSGE.DSGE_DATE_FORMAT)
    @test get_setting(m, :cond_full_names) ==  [:obs_gdp, :obs_corepce, :obs_spread,
                                                            :obs_nominalrate, :obs_longrate]
    @test get_setting(m, :cond_semi_names) == [:obs_spread, :obs_nominalrate, :obs_longrate]
    @test get_setting(m, :use_population_forecast) == false
    @test get_setting(m, :population_mnemonic).value == :CNP16OV__FRED
    @test get_setting(m, :hpfilter_population) == true
    @test get_setting(m, :rate_expectations_source) == :ois

    # Dates
    @test get_setting(m, :date_presample_start) == DSGE.quartertodate("1959-Q3")
    @test get_setting(m, :date_mainsample_start) == DSGE.quartertodate("1960-Q1")
    @test get_setting(m, :date_zlb_start) == DSGE.quartertodate("2008-Q4")
    @test get_setting(m, :date_forecast_start) == Dates.lastdayofquarter(Dates.today())
    @test get_setting(m, :date_conditional_end) == Dates.lastdayofquarter(Dates.today())

    # Anticipated shocks
    @test get_setting(m, :n_mon_anticipated_shocks) == 0
    @test get_setting(m, :n_mon_anticipated_shocks_padding) == 20

    # General computation
    @test get_setting(m, :use_parallel_workers) == true

    # Optimization
    @test get_setting(m, :reoptimize) == true
    @test get_setting(m, :optimization_method) == :csminwel
    @test get_setting(m, :optimization_iterations) == 100
    @test get_setting(m, :optimization_step_size) == .01
    @test get_setting(m, :simulated_annealing_temperature) == Optim.log_temperature
    @test get_setting(m, :simulated_annealing_block_proportion) == .3
    @test get_setting(m, :optimization_ftol) == 1e-10
    @test get_setting(m, :optimization_xtol) == 1e-10
    @test get_setting(m, :optimization_gtol) == 1e-10
    @test get_setting(m, :combined_optimizer_max_cycles) == 4
    @test get_setting(m, :optimization_attempts) == 4

    # Estimation
    @test get_setting(m, :sampling_method) == :MH
    @test get_setting(m, :calculate_hessian) == true
    @test get_setting(m, :n_hessian_test_params) == typemax(Int)
    @test get_setting(m, :use_chand_recursion) == false

    # Metropolis-Hastings
    @test get_setting(m, :n_mh_simulations) == 5000
    @test get_setting(m, :n_mh_blocks) == 22
    @test get_setting(m, :n_mh_param_blocks) == 1
    @test get_setting(m, :mh_adaptive_accept) == false
    @test get_setting(m, :mh_target_accept) == .25
    @test get_setting(m, :mh_c) == 0.5
    @test get_setting(m, :mh_α) == 1.0
    @test get_setting(m, :n_mh_burn) == 2
    @test get_setting(m, :mh_thin) == 5
    @test get_setting(m, :mh_cc) == 0.09
    @test get_setting(m, :mh_cc0) == 0.01

    # Forecast
    @test get_setting(m, :forecast_block_size) == 5000
    @test get_setting(m, :forecast_start_block) == 1
    @test get_setting(m, :forecast_input_file_overrides) == Dict{Symbol,String}()
    @test get_setting(m, :forecast_jstep) == 5
    @test get_setting(m, :forecast_uncertainty_override).hasvalue == false
    @test get_setting(m, :forecast_smoother) == :durbin_koopman
    @test get_setting(m, :forecast_horizons) == 60
    @test get_setting(m, :forecast_tdist_shocks) == false
    @test get_setting(m, :forecast_tdist_df_val) == 15
    @test get_setting(m, :forecast_zlb_value) == 0.13 / 4
    @test get_setting(m, :shockdec_startdate).hasvalue == false
    @test get_setting(m, :shockdec_enddate).hasvalue == false
    @test get_setting(m, :impulse_response_horizons) == 40
    @test get_setting(m, :compute_shockdec_bands) == false

    # Sequential Monte Carlo
    @test get_setting(m, :n_particles) == 2000
    @test get_setting(m, :n_Φ) == 300
    @test get_setting(m, :λ) == 2.1
    @test get_setting(m, :n_smc_blocks) == 1
    @test get_setting(m, :step_size_smc) == .5
    @test get_setting(m, :n_mh_steps_smc) == 1
    @test get_setting(m, :target_accept) == 0.25
    @test get_setting(m, :resampler_smc) == :multinomial
    @test get_setting(m, :mixture_proportion) == 1.


    # Adaptive ϕ Schedule
    @test get_setting(m, :use_fixed_schedule) == true
    @test get_setting(m, :tempering_target) == 0.95
    @test get_setting(m, :resampling_threshold) == 0.5

    # Temporary setting to save different output files
    @test get_setting(m, :adaptive_tempering_target_smc) == 0.97
    @test get_setting(m, :smc_iteration) == 1
    @test get_setting(m, :previous_data_vintage) == Dates.format(now(), DSGE.DSGE_DATE_FORMAT)

    # Alternative policy
    @test !haskey(DSGE.get_settings(m), :alternative_policy)
end

m = AnSchorfheide(testing = true)
@testset "Check default test settings" begin
    # I/O
    @test get_setting(m, :dataroot) == normpath(joinpath(dirname(@__FILE__), "..", "test", "reference", "input_data"))

    # General
    @test get_setting(m, :data_vintage) == "REF"
    @test get_setting(m, :cond_vintage) == "REF"
    @test get_setting(m, :use_parallel_workers) == false
    @test get_setting(m, :n_hessian_test_params) == 3

    # Metropolis-Hastings
    @test get_setting(m, :n_mh_simulations) == 100
    @test get_setting(m, :n_mh_blocks) == 1
    @test get_setting(m, :n_mh_burn) == 0
    @test get_setting(m, :mh_thin) == 1

    # Forecast
    @test get_setting(m, :date_forecast_start) == DSGE.quartertodate("2015-Q4")
    @test get_setting(m, :forecast_horizons) == 2
    @test get_setting(m, :forecast_jstep) == 1
    @test get_setting(m, :impulse_response_horizons) == 2
end
