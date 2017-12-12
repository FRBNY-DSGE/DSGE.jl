"""
```
default_settings!(m::AbstractModel)
```

Default Settings are constructed, initialized and added to `m.settings`.
"""
function default_settings!(m::AbstractModel)

    settings = m.settings

    # I/O File locations
    saveroot = normpath(joinpath(dirname(@__FILE__), "..","save"))
    datapath = normpath(joinpath(dirname(@__FILE__), "..","save","input_data"))

    settings[:saveroot] = Setting(:saveroot, saveroot, "Root of data directory structure")
    settings[:dataroot] = Setting(:dataroot, datapath, "Input data directory path")

    # Data settings for released and conditional data. Default behavior is to set vintage
    # of data to today's date.
    vint = Dates.format(now(), DSGE_DATE_FORMAT)
    settings[:data_vintage] = Setting(:data_vintage, vint, true, "vint",
        "Data vintage")
    settings[:data_id] = Setting(:data_id, 3,
        "Dataset identifier")
    settings[:cond_vintage] = Setting(:cond_vintage, vint,
        "Conditional data vintage")
    settings[:cond_id] = Setting(:cond_id, 2,
        "Conditional dataset identifier")
    settings[:cond_full_names] = Setting(:cond_full_names, [:obs_gdp, :obs_corepce, :obs_spread, :obs_nominalrate, :obs_longrate],
        "Observables used in conditional forecasts")
    settings[:cond_semi_names] = Setting(:cond_semi_names, [:obs_spread, :obs_nominalrate, :obs_longrate],
        "Observables used in semiconditional forecasts")
    settings[:use_population_forecast] = Setting(:use_population_forecast, false,
        "Whether to use population forecasts as data")
    settings[:population_mnemonic] = Setting(:population_mnemonic, Nullable(:CNP16OV__FRED),
        "Mnemonic of FRED data series for computing per-capita values (a Nullable{Symbol})")
    settings[:hpfilter_population] = Setting(:hpfilter_population, true,
        "Whether to HP filter combined population and forecast")

    # Dates
    settings[:date_presample_start] = Setting(:date_presample_start, quartertodate("1959-Q3"),
        "Start date of pre-sample")
    settings[:date_mainsample_start] = Setting(:date_mainsample_start, quartertodate("1960-Q1"),
        "Start date of main sample")
    settings[:date_zlb_start] = Setting(:date_zlb_start, quartertodate("2008-Q4"),
        "Start date of zero lower bound regime")
    settings[:date_forecast_start] = Setting(:date_forecast_start, Dates.lastdayofquarter(Dates.today()),
        "Start date of forecast period")
    settings[:date_conditional_end] = Setting(:date_conditional_end, Dates.lastdayofquarter(Dates.today()),
        "End date of conditional data period")

    # Anticipated shocks
    settings[:n_anticipated_shocks] = Setting(:n_anticipated_shocks, 0,
        "Number of anticipated policy shocks")
    settings[:n_anticipated_shocks_padding] = Setting(:n_anticipated_shocks_padding, 20,
        "Padding for anticipated policy shocks")

    # General computation
    settings[:use_parallel_workers] = Setting(:use_parallel_workers, true,
        "Use available parallel workers in computations")

    # Optimization
    settings[:reoptimize] = Setting(:reoptimize, true,
        "Optimize the posterior mode. If false, reads in mode from a file.")
    settings[:optimization_method] = Setting(:optimization_method, :csminwel,
        "Method for finding the posterior mode")
    settings[:optimization_iterations] = Setting(:optimization_iterations, 100,
        "Number of iterations the optimizer should run for")
    settings[:optimization_step_size] = Setting(:optimization_step_size, 0.01,
        "Step size scaling factor for optimization")
	settings[:simulated_annealing_temperature] = Setting(:simulated_annealing_temperature, Optim.log_temperature,
        "Temperature function for simulated annealing")
   settings[:simulated_annealing_block_proportion] = Setting(:simulated_annealing_block_proportion, 0.3,
        "Fraction of parameters to vary for each proposed move in simulated annealing")
   settings[:optimization_ftol] = Setting(:optimization_ftol, 1e-10,
        "Relative function difference threshold for optimization")
    settings[:optimization_xtol] = Setting(:optimization_xtol, 1e-10,
        "Relative input vector difference threshold for optimization")
    settings[:optimization_gtol] = Setting(:optimization_gtol, 1e-10,
        "Relative gradient difference threshold for optimization")
    settings[:combined_optimizer_max_cycles] = Setting(:combined_optimizer_max_cycles, 4,
        "Total number of cycles to use in the combined optimization routine")
    settings[:optimization_attempts] = Setting(:optimization_attempts, 4,
        "Number of times to attempt optimization in estimate()")

    # Estimation
    settings[:sampling_method] = Setting(:sampling_method, :MH, "The method used to generate samples from the posterior distribution of the parameters.")
    settings[:calculate_hessian] = Setting(:calculate_hessian, true,
        "Calculate the hessian at the mode")
    settings[:n_hessian_test_params] = Setting(:n_hessian_test_params, typemax(Int),
        "Max number of free params for which to calculate Hessian")

    # Metropolis-Hastings
    settings[:n_mh_simulations] = Setting(:n_mh_simulations, 5000,
        "Number of draws saved (after thinning) per block in Metropolis-Hastings")
    settings[:n_mh_blocks] = Setting(:n_mh_blocks, 22,
        "Number of blocks for Metropolis-Hastings")
    settings[:n_mh_burn] = Setting(:n_mh_burn, 2,
        "Number of blocks to use as burn-in in Metropolis-Hastings")
    settings[:mh_thin] = Setting(:mh_thin, 5,
        "Metropolis-Hastings thinning step")

    # Forecast
    settings[:forecast_block_size] = Setting(:forecast_block_size, 5000,
        "Number of draws in each forecast block (before thinning by forecast_jstep)")
    settings[:forecast_start_block] = Setting(:forecast_start_block, 1,
        "Block at which to resume forecasting")
    settings[:forecast_input_file_overrides] = Setting(:forecast_input_file_overrides,
        Dict{Symbol, String}())
    settings[:forecast_jstep] = Setting(:forecast_jstep, 5,
        "Forecast thinning step (in addition to MH thinning step")
    settings[:forecast_uncertainty_override] = Setting(:forecast_uncertainty_override, Nullable{Bool}(),
        "If non-null, overrides default drawing states/shocks behavior in smoother and forecast")
    settings[:forecast_smoother] = Setting(:forecast_smoother, :durbin_koopman,
        "Choice of smoother to use during forecasting. Can be :hamilton, :koopman, :carter_kohn, or :durbin_koopman")
    settings[:forecast_horizons] = Setting(:forecast_horizons, 60,
        "Number of periods to forecast ahead")
    settings[:forecast_tdist_shocks] = Setting(:forecast_tdist_shocks, false,
        "Draw Students-t distributed shocks in forecast")
    settings[:forecast_tdist_df_val] = Setting(:forecast_tdist_df_val, 15,
        "Students-t degrees of freedom fixed value")
    settings[:forecast_zlb_value] = Setting(:forecast_zlb_value, 0.13/4,
        "Value of the zero lower bound in forecast periods, if we choose to enforce it")
    settings[:shockdec_startdate] = Setting(:shockdec_startdate, Nullable{Date}(),
        "Date of start of shock decomposition output period. If null, then shockdec starts at date_mainsample_start")
    settings[:shockdec_enddate] = Setting(:shockdec_enddate, Nullable{Date}(),
        "Date of end of shock decomposition output period. If null, then shockdec ends at date_forecast_end")
    settings[:impulse_response_horizons] = Setting(:impulse_response_horizons, 40,
        "Number of periods for which to calculate an impulse response")
    settings[:compute_shockdec_bands] = Setting(:compute_shockdec_bands, false, "Whether or not to compute bands for shock decomposition. Setting to false saves signficant storage space.")

	# Sequential Monte Carlo
    settings[:n_particles] = Setting(:n_particles, 2000, "Number of particles for use in SMC")
    settings[:n_Φ] = Setting(:n_Φ, 300, "Number of stages in the tempering schedule")
  	settings[:λ] = Setting(:λ, 2.1, "The 'bending coefficient' λ in Φ(n) = (n/N(Φ))^λ")
    settings[:n_smc_blocks] = Setting(:n_smc_blocks, 1, "The number of parameter blocks in SMC")
    settings[:step_size_smc] = Setting(:step_size_smc, .5, "The scaling factor for the covariance of the particles. Controls size of steps in mutation step")
    settings[:n_MH_steps_smc] = Setting(:n_MH_steps_smc, 5, "Number of Metropolis Hastings steps to attempt during the mutation step.")
    settings[:init_accept] = Setting(:init_accept, .25, "The initial average acceptance rate for new particles during mutation")
    settings[:target_accept] = Setting(:target_accept, .25, "The initial target acceptance rate for new particles during mutation")
    settings[:resampler_smc] = Setting(:resampler_smc, :multinomial, "Which resampling method to use in SMC")
    settings[:initial_draw_source] = Setting(:initial_draw_source, :prior, "How to draw the initial population of particles in SMC")

    # Endogenous ϕ Schedule
    settings[:use_fixed_schedule] = Setting(:use_fixed_schedule, true, true, "fix", "Boolean indicating whether or not to use a fixed tempering (ϕ) schedule")
    settings[:tempering_target] = Setting(:tempering_target, 0.95, "The coefficient of the sample size metric to be targeted when solving for an endogenous ϕ")
    settings[:resampling_threshold] = Setting(:resampling_threshold, 0.5, "The threshold such that the particles will be resampled when the population drops below threshold * N")
    # temporary setting to save different output files
    settings[:smc_iteration] = Setting(:smc_iteration, 1, true, "iter", "The iteration index for the number of times smc has been run on the same data vintage. Primarily for numerical accuracy/testing purposes.")

    # Time tempering (to be changed)
    settings[:updated_data_vintage] = Setting(:updated_data_vintage, vint, "Updated data vintage for the combination tempered data/likelihood
                                              tempering for the purposes of updating an old cloud vintage with new data.")
    settings[:date_updatedsample_start] = Setting(:date_updatedsample_start, quartertodate("1960-Q1"),
                                                 "Start date of revision or update of sample")
    settings[:date_updatedsample_end] = Setting(:date_updatedsample_end, quartertodate("1960-Q2"),
                                                 "End date of revision or update of sample")

    # Alternative policy
    baseline_policy = AltPolicy(:historical, eqcond, solve, forecast_init = identity)
    settings[:alternative_policy] = Setting(:alternative_policy, baseline_policy)

    return settings
end

"""
```
default_test_settings!(m::AbstractModel)
```

The following Settings are constructed, initialized and added to
`m.test_settings`. Their purposes are identical to those in
`m.settings`, but these values are used to test DSGE.jl.

### I/O Locations and identifiers
- `saveroot::Setting{String}`: A temporary directory in /tmp/
- `dataroot::Setting{String}`: dsgeroot/test/reference/
- `data_vintage::Setting{String}`: \"_REF\"

### Metropolis-Hastings
- `n_mh_simulations::Setting{Int}`: 100
- `n_mh_blocks::Setting{Int}`: 1
- `n_mh_burn::Setting{Int}`: 0
- `mh_thin::Setting{Int}`: 1
"""
function default_test_settings!(m::AbstractModel)

    test = m.test_settings

    # I/O
    dataroot = normpath(joinpath(dirname(@__FILE__), "..", "test", "reference", "input_data"))
    saveroot = mktempdir()

    # General
    test[:saveroot] = Setting(:saveroot, saveroot,
        "Where to write files when in test mode")
    test[:dataroot] = Setting(:dataroot, dataroot,
        "Location of input files when in test mode" )
    test[:data_vintage] = Setting(:data_vintage, "REF", true, "vint",
        "Reference data identifier")
    test[:cond_vintage] = Setting(:cond_vintage, "REF",
        "Vintage of conditional data")
    test[:use_parallel_workers] = Setting(:use_parallel_workers, false, false, "parw",
        "Use available parallel workers in computations")
    test[:n_hessian_test_params] = Setting(:n_hessian_test_params, 3, false, "mhfp",
        "Max number of free params for which to calculate Hessian")

    # Metropolis-Hastings
    test[:n_mh_simulations] = Setting(:n_mh_simulations, 100, false, "nsim",
        "Number of parameter draws per block for testing Metropolis-Hastings")
    test[:n_mh_blocks] = Setting(:n_mh_blocks, 1, false, "nblc",
        "Number of blocks to draw parameters for testing Metropolis-Hastings")
    test[:n_mh_burn] = Setting(:n_mh_burn, 0, false, "nbrn",
        "Number of burn-in blocks for testing Metropolis-Hastings")
    test[:mh_thin] = Setting(:mh_thin, 1, false, "thin",
        "Thinning step for testing Metropolis-Hastings")

    # Forecast
    test[:date_forecast_start] = Setting(:date_forecast_start, quartertodate("2015-Q4"),
        "Start date of forecast period")
    test[:forecast_horizons] = Setting(:forecast_horizons, 2,
        "Number of periods to forecast ahead")
    test[:forecast_jstep] = Setting(:forecast_jstep, 1,
        "Forecast thinning step (in addition to MH thinning step")
    test[:impulse_response_horizons] = Setting(:impulse_response_horizons, 2,
        "Number of periods for which to calculate an impulse response")

    return test
end
