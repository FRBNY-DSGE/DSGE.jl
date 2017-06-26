isdefined(Base, :__precompile__) && __precompile__()

module DSGE
    using Base.Dates, DataFrames, Distributions, FredData, HDF5, JLD, Optim, StateSpaceRoutines
    using DataStructures: SortedDict, insert!, ForwardOrdering, OrderedDict
    using QuantEcon: solve_discrete_lyapunov
    using Roots: fzero, ConvergenceFailed
    using Base.Test
    import Calculus
    import Optim: optimize, Optimizer

    export

        # distributions_ext.jl
        BetaAlt, GammaAlt, DegenerateMvNormal, DegenerateDiagMvTDist,

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
        n_parameters_free, n_pseudoobservables, get_key,
        inds_states_no_ant, inds_shocks_no_ant, inds_obs_no_ant,
        spec, subspec, saveroot, dataroot,
        data_vintage, cond_vintage, cond_id, cond_full_names, cond_semi_names, use_population_forecast,
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

        # parameters.jl
        parameter, Transform, NullablePrior, AbstractParameter,
        Parameter, ParameterVector, ScaledParameter,
        UnscaledParameter, SteadyStateParameter, transform_to_real_line, transform_to_model_space,
        update, update!, transform_to_model_space, transform_to_real_line, Interval, ParamBoundsError,

        # observables.jl
        Observable, PseudoObservable, PseudoObservableMapping, check_mnemonics,

        # statespace.jl
        Measurement, Transition, System, compute_system,

        # estimate/
        simulated_annealing, combined_optimizer, LBFGS_wrapper,
        filter, likelihood, posterior, posterior!,
        optimize!, csminwel, hessian!, estimate, proposal_distribution,

        metropolis_hastings, compute_parameter_covariance, compute_moments,
        find_density_bands, prior, mutation, multinomial_resampling, systematic_resampling, mutation, smc, nearestSPD, density, correction, ineff_func,

        # forecast/
        load_draws, forecast_one,
        smooth, forecast, shock_decompositions, deterministic_trends, trends, impulse_responses,
        compute_system, add_requisite_output_vars, n_forecast_draws,
        get_forecast_input_file, get_forecast_output_files, get_forecast_filename,
        read_forecast_output,

        # models/
        init_parameters!, steadystate!, init_observable_mappings!,
        Model990, Model1002, Model1010, SmetsWouters, AnSchorfheide, eqcond, measurement, pseudo_measurement,

        # solve/
        gensys, solve,

        # data/
        load_data, load_data_levels, load_cond_data_levels, load_fred_data,
        transform_data, save_data, get_data_filename,
        df_to_matrix, hpfilter, difflog, quartertodate, percapita, nominal_to_real,
        oneqtrpctchange, annualtoquarter, quartertoannual, quartertoannualpercent,
        loggrowthtopct_annualized_percapita, loggrowthtopct_annualized, logleveltopct_annualized_percapita,
        logleveltopct_annualized,
        parse_data_series, collect_data_transforms, reverse_transform,
        subtract_quarters, iterate_quarters,

        # analysis/
        find_density_bands, moment_tables, means_bands, means_bands_all, compute_means_bands, MeansBands,
        meansbands_matrix_all, meansbands_matrix, read_mb,
        get_meansbands_input_files, get_meansbands_output_files, get_product, get_class,
        which_density_bands, write_meansbands_tables, prepare_meansbands_tables_timeseries,
        prepare_meansbands_tables_shockdec, write_meansbands_tables_all,

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
    include("estimate/LBFGS.jl")
    include("estimate/estimate.jl")
    include("estimate/mutation_RWMH.jl")
    include("estimate/systematic_resampling.jl")
    include("estimate/smc.jl")
    include("estimate/nearestSPD.jl")
    include("estimate/multinomial_resampling.jl")
    include("estimate/mutation.jl")
    include("estimate/correction.jl")
    include("estimate/ineff_func.jl")
    include("estimate/density.jl")

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

    include("models/financial_frictions.jl")

    include("models/m990/m990.jl")
    include("models/m990/subspecs.jl")
    include("models/m990/eqcond.jl")
    include("models/m990/observables.jl")
    include("models/m990/measurement.jl")
    include("models/m990/pseudo_measurement.jl")
    include("models/m990/augment_states.jl")

    include("models/m1002/m1002.jl")
    include("models/m1002/subspecs.jl")
    include("models/m1002/eqcond.jl")
    include("models/m1002/observables.jl")
    include("models/m1002/measurement.jl")
    include("models/m1002/pseudo_measurement.jl")
    include("models/m1002/augment_states.jl")

    include("models/m1010/m1010.jl")
    include("models/m1010/subspecs.jl")
    include("models/m1010/eqcond.jl")
    include("models/m1010/observables.jl")
    include("models/m1010/measurement.jl")
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
    include("models/an_schorfheide/pseudo_measurement.jl")
    include("models/an_schorfheide/augment_states.jl")
end
