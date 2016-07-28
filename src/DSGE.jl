isdefined(Base, :__precompile__) && __precompile__()

module DSGE
    using Distributions, Roots.fzero, HDF5, JLD
    using DataStructures: SortedDict, insert!, ForwardOrdering, OrderedDict
    using FredData, DataFrames, Base.Dates
    using QuantEcon: solve_discrete_lyapunov
    import Calculus
    import Optim
    using Optim: OptimizationTrace, OptimizationState, MultivariateOptimizationResults
    import NaNMath

    export

        # distributions_ext.jl
        BetaAlt, GammaAlt, DegenerateMvNormal,

        # settings.jl
        Setting, get_setting, default_settings!, default_test_settings!,

        # abstractdsgemodel.jl
        AbstractModel,
        n_anticipated_shocks, n_anticipated_shocks_padding,
        date_presample_start, date_prezlb_start, date_zlb_start,
            date_presample_end, date_prezlb_end, date_zlb_end,
            index_presample_start, index_prezlb_start, index_zlb_start,
            n_presample_periods, n_prezlb_periods, n_zlb_periods,
        n_states, n_states_augmented, n_shocks_exogenous, n_shocks_expectational,
            n_equilibrium_conditions, n_observables, n_parameters, n_parameters_steady_state,
            n_parameters_free,
        spec, subspec, saveroot, dataroot,
        data_vintage, cond_vintage, cond_id, use_population_forecast,
        use_parallel_workers,
        reoptimize, calculate_hessian, hessian_path, n_hessian_test_params,
        n_mh_blocks, n_mh_simulations, n_mh_burn, mh_thin, n_draws,
        date_forecast_start, forecast_tdist_df_val, forecast_tdist_shocks, forecast_kill_shocks,
            forecast_smoother, forecast_horizons,
        load_parameters_from_file, specify_mode!, specify_hessian,
        logpath, workpath, rawpath, tablespath, figurespath, inpath,
        transform_to_model_space!, transform_to_real_line!,

        # parameters.jl
        parameter, Transform, NullablePrior, AbstractParameter,
        Parameter, ParameterVector, ScaledParameter,
        UnscaledParameter, SteadyStateParameter, transform_to_real_line, transform_to_model_space,
        update, update!, transform_to_model_space, transform_to_real_line, Interval, ParamBoundsError,

        # statespace/
        Measurement, Transition, System, compute_system,

        # estimate/
        kalman_filter, kalman_filter_2part, likelihood, posterior, posterior!,
        optimize!, csminwel, hessian!, estimate, proposal_distribution,
        metropolis_hastings, compute_parameter_covariance, compute_moments,
        find_density_bands, prior,

        # forecast/
        kalman_smoother, disturbance_smoother, forecast_all, forecast_one, computeForecast,
        durbin_koopman_smoother,

        # models/
        steadystate!, Model990, SmetsWouters, eqcond, measurement,

        # solve/
        gensys, solve,

        # data/
        load_fred_data, load_data, transform_data, save_data_matrix!, hpfilter, difflog,
        quartertodate, df_to_matrix

    const VERBOSITY = Dict(:none => 0, :low => 1, :high => 2)
    const DSGE_DATE_FORMAT = "yymmdd"

    include("parameters.jl")
    include("distributions_ext.jl")
    include("abstractdsgemodel.jl")
    include("settings.jl")
    include("defaults.jl")
    include("statespace.jl")
    include("util.jl")

    include("data/load_data.jl")
    include("data/fred_data.jl")
    include("data/transformations.jl")
    include("data/transform_data.jl")
    include("data/load_cond_data.jl")
    include("data/util.jl")

    include("solve/gensys.jl")
    include("solve/solve.jl")

    include("estimate/kalman.jl")
    include("estimate/posterior.jl")
    include("estimate/optimize.jl")
    include("estimate/csminwel.jl")
    include("estimate/hessian.jl")
    include("estimate/hessizero.jl")
    include("estimate/estimate.jl")
    include("estimate/moments.jl")

    include("forecast/smoothers.jl")
    include("forecast/compute_system.jl")

    include("forecast/filter.jl")
    include("forecast/forecast_all.jl")
    include("forecast/forecast.jl")
    include("forecast/computeForecast.jl")

    include("models/m990/m990.jl")
    include("models/m990/subspecs.jl")
    include("models/m990/eqcond.jl")
    include("models/m990/measurement.jl")
    include("models/m990/augment_states.jl")

    include("models/smets_wouters/smets_wouters.jl")
    include("models/smets_wouters/subspecs.jl")
    include("models/smets_wouters/eqcond.jl")
    include("models/smets_wouters/measurement.jl")
    include("models/smets_wouters/augment_states.jl")

    include("data/load_data.jl")
    include("data/fred_data.jl")
    include("data/transform_data.jl")
    include("data/transformations.jl")
    include("data/util.jl")

end
