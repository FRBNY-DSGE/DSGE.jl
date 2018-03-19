isdefined(Base, :__precompile__) && __precompile__()

module DSGE
    using Base.Dates, Base.Test
    using DataFrames, Distributions, FredData, HDF5, JLD, Optim, Plots, RecipesBase, StateSpaceRoutines, StatPlots
    using DataStructures: SortedDict, insert!, ForwardOrdering, OrderedDict
    using QuantEcon: solve_discrete_lyapunov
    using Roots: fzero, ConvergenceFailed
    using StatsBase: sample
    import Calculus
    import Optim: optimize, SecondOrderOptimizer, MultivariateOptimizationResults

    export

        # distributions_ext.jl
        BetaAlt, GammaAlt, DegenerateMvNormal, DegenerateDiagMvTDist, MatrixNormal,

        # settings.jl
        Setting, get_setting,

        # defaults.jl
        default_settings!, default_test_settings!,

        # abstractdsgemodel.jl
        AbstractModel, description,
        n_anticipated_shocks, n_anticipated_shocks_padding,
        date_presample_start, date_mainsample_start, date_zlb_start,
        date_presample_end, date_prezlb_end, date_mainsample_end, date_conditional_end,
        index_presample_start, index_mainsample_start, index_zlb_start, index_forecast_start,
        index_shockdec_start,
        n_presample_periods, n_prezlb_periods, n_zlb_periods, n_mainsample_periods,
        n_conditional_periods,
        inds_presample_periods, inds_prezlb_periods, inds_zlb_periods, inds_mainsample_periods,
        n_states, n_states_augmented, n_shocks_exogenous, n_shocks_expectational,
        n_equilibrium_conditions, n_observables, n_parameters, n_parameters_steady_state,
        n_parameters_free, n_pseudo_observables, get_key,
        inds_states_no_ant, inds_shocks_no_ant, inds_obs_no_ant,
        spec, subspec, saveroot, dataroot,
        data_vintage, data_id, cond_vintage, cond_id, cond_full_names, cond_semi_names, use_population_forecast,
        use_parallel_workers,
        reoptimize, calculate_hessian, hessian_path, n_hessian_test_params,
        n_mh_blocks, n_mh_simulations, n_mh_burn, mh_thin,
        date_forecast_start, date_forecast_end,
        forecast_block_size, forecast_start_block,
        forecast_input_file_overrides, forecast_uncertainty_override,
        forecast_smoother, forecast_horizons,
        forecast_zlb_value, forecast_tdist_shocks, forecast_tdist_df_val,
        shockdec_startdate, date_shockdec_end,
        n_shockdec_periods, impulse_response_horizons,
        load_parameters_from_file, specify_mode!, specify_hessian,
        logpath, workpath, rawpath, tablespath, figurespath, inpath,
        transform_to_model_space!, transform_to_real_line!,
        ShockGroup, alternative_policy,

        # parameters.jl
        parameter, Transform, NullablePrior, AbstractParameter,
        Parameter, ParameterVector, ScaledParameter,
        UnscaledParameter, SteadyStateParameter, transform_to_real_line, transform_to_model_space,
        update, update!, transform_to_model_space, transform_to_real_line, Interval, ParamBoundsError,

        # observables.jl
        Observable, PseudoObservable, check_mnemonics,

        # statespace.jl
        Transition, Measurement, PseudoMeasurement, System, compute_system,

        # data/
        load_data, load_data_levels, load_cond_data_levels, load_fred_data,
        transform_data, save_data, get_data_filename,
        df_to_matrix, hpfilter, difflog, quartertodate, percapita, nominal_to_real,
        oneqtrpctchange, annualtoquarter, quartertoannual, quartertoannualpercent,
        loggrowthtopct_percapita, loggrowthtopct, logleveltopct_annualized,
        loggrowthtopct_annualized_percapita, loggrowthtopct_annualized, logleveltopct_annualized_percapita,
        logleveltopct_annualized_approx, loggrowthtopct_4q_approx, logleveltopct_4q_approx,
        parse_data_series, collect_data_transforms, reverse_transform,
        subtract_quarters, iterate_quarters,

        # solve/
        gensys, solve,

        # estimate/
        simulated_annealing, combined_optimizer, lbfgs,
        filter, likelihood, posterior, posterior!,
        optimize!, csminwel, hessian!, estimate, proposal_distribution,
        metropolis_hastings, compute_parameter_covariance,
        prior, get_estimation_output_files,

        # forecast/
        load_draws, forecast_one,
        smooth, forecast, shock_decompositions, deterministic_trends, trends, impulse_responses,
        compute_system, compute_system_function, add_requisite_output_vars, n_forecast_draws,
        get_forecast_input_file, get_forecast_output_files, get_forecast_filename,
        read_forecast_output,

        # analysis/
        find_density_bands, moment_tables, compute_meansbands, MeansBands,
        meansbands_to_matrix, read_mb, read_bdd_and_unbdd_mb,
        get_meansbands_input_file, get_meansbands_output_file, get_product, get_class,
        which_density_bands,
        prepare_meansbands_tables_timeseries, prepare_means_tables_shockdec, prepare_meansbands_table_irf,
        write_meansbands_tables_timeseries, write_means_tables_shockdec, prepare_meansbands_table_irf,
        write_meansbands_tables_all,
        construct_fcast_and_hist_dfs,
        df_to_table,

        # altpolicy/
        AltPolicy, taylor93, taylor99,

        # scenarios/
        AbstractScenario, SingleScenario, Scenario, SwitchingScenario, ScenarioAggregate,
        n_targets, n_instruments, n_target_horizons, targets_to_data,
        compute_scenario_system, filter_shocks!, forecast_scenario, simulate_switching, scenario_means_bands,
        get_scenario_input_file, n_scenario_draws, get_scenario_filename, get_scenario_output_files,
        read_scenario_output, get_scenario_mb_input_file, get_scenario_mb_output_file, read_scenario_mb,

        # plot/
        plot_prior_posterior, plot_impulse_response, plot_history_and_forecast, hair_plot,
        plot_forecast_comparison, plot_shock_decomposition, plot_altpolicies, plot_scenario,

        # models/
        init_parameters!, steadystate!, init_observable_mappings!, init_pseudo_observable_mappings!,
        Model990, Model1002, Model1010, SmetsWouters, AnSchorfheide, eqcond, measurement, pseudo_measurement,
        shock_groupings,

        # util
        @test_matrix_approx_eq, @test_matrix_approx_eq_eps

    const VERBOSITY = Dict(:none => 0, :low => 1, :high => 2)
    const DSGE_DATE_FORMAT = "yymmdd"
    const DSGE_DATASERIES_DELIM = "__"
    const DSGE_SHOCKDEC_DELIM = "__"

    include("parameters.jl")
    include("distributions_ext.jl")
    include("abstractdsgemodel.jl")
    include("settings.jl")
    include("defaults.jl")
    include("observables.jl")
    include("statespace.jl")
    include("util.jl")

    include("data/load_data.jl")
    include("data/fred_data.jl")
    include("data/transformations.jl")
    include("data/transform_data.jl")
    include("data/reverse_transform.jl")
    include("data/util.jl")

    include("solve/gensys.jl")
    include("solve/solve.jl")

    include("estimate/kalman.jl")
    include("estimate/filter.jl")
    include("estimate/posterior.jl")
    include("estimate/optimize.jl")
    include("estimate/csminwel.jl")
    include("estimate/hessian.jl")
    include("estimate/hessizero.jl")
    include("estimate/simulated_annealing.jl")
    include("estimate/combined_optimizer.jl")
    include("estimate/lbfgs.jl")
    include("estimate/nelder_mead.jl")
    include("estimate/estimate.jl")

    include("forecast/util.jl")
    include("forecast/io.jl")
    include("forecast/smooth.jl")
    include("forecast/forecast.jl")
    include("forecast/shock_decompositions.jl")
    include("forecast/impulse_responses.jl")
    include("forecast/drivers.jl")

    include("analysis/moments.jl")
    include("analysis/meansbands.jl")
    include("analysis/compute_meansbands.jl")
    include("analysis/meansbands_to_matrix.jl")
    include("analysis/io.jl")
    include("analysis/util.jl")
<<<<<<< HEAD
=======
    include("analysis/forecast_decomposition.jl")
    include("analysis/df_to_table.jl")
>>>>>>> 28fef3f... Create functions in analysis/df_to_table for producing LaTeX tables from dfs

    include("altpolicy/altpolicy.jl")
    include("altpolicy/taylor93.jl")
    include("altpolicy/taylor99.jl")

    include("scenarios/scenario.jl")
    include("scenarios/io.jl")
    include("scenarios/forecast.jl")
    include("scenarios/switching.jl")
    include("scenarios/transform.jl")

    include("plot/util.jl")
    include("plot/plot_prior_posterior.jl")
    include("plot/plot_impulse_response.jl")
    include("plot/plot_history_and_forecast.jl")
    include("plot/hair_plot.jl")
    include("plot/plot_forecast_comparison.jl")
    include("plot/plot_shock_decomposition.jl")
    include("plot/plot_altpolicies.jl")
    include("plot/plot_scenario.jl")

    include("models/financial_frictions.jl")

    include("models/m990/m990.jl")
    include("models/m990/subspecs.jl")
    include("models/m990/eqcond.jl")
    include("models/m990/observables.jl")
    include("models/m990/measurement.jl")
    include("models/m990/pseudo_observables.jl")
    include("models/m990/pseudo_measurement.jl")
    include("models/m990/augment_states.jl")

    include("models/m1002/m1002.jl")
    include("models/m1002/subspecs.jl")
    include("models/m1002/eqcond.jl")
    include("models/m1002/observables.jl")
    include("models/m1002/measurement.jl")
    include("models/m1002/pseudo_observables.jl")
    include("models/m1002/pseudo_measurement.jl")
    include("models/m1002/augment_states.jl")

    include("models/m1010/m1010.jl")
    include("models/m1010/subspecs.jl")
    include("models/m1010/eqcond.jl")
    include("models/m1010/observables.jl")
    include("models/m1010/measurement.jl")
    include("models/m1010/pseudo_observables.jl")
    include("models/m1010/pseudo_measurement.jl")
    include("models/m1010/augment_states.jl")

    include("models/smets_wouters/smets_wouters.jl")
    include("models/smets_wouters/subspecs.jl")
    include("models/smets_wouters/eqcond.jl")
    include("models/smets_wouters/observables.jl")
    include("models/smets_wouters/measurement.jl")
    include("models/smets_wouters/augment_states.jl")

    include("models/an_schorfheide/an_schorfheide.jl")
    include("models/an_schorfheide/subspecs.jl")
    include("models/an_schorfheide/eqcond.jl")
    include("models/an_schorfheide/observables.jl")
    include("models/an_schorfheide/measurement.jl")
    include("models/an_schorfheide/pseudo_observables.jl")
    include("models/an_schorfheide/pseudo_measurement.jl")
    include("models/an_schorfheide/augment_states.jl")
end
