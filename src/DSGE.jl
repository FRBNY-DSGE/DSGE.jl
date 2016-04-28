isdefined(Base, :__precompile__) && __precompile__()

module DSGE
    using Distributions, Roots.fzero, HDF5
    using DataStructures: SortedDict, insert!, ForwardOrdering, OrderedDict
    using FredData, DataFrames, Base.Dates
    using QuantEcon: solve_discrete_lyapunov
    import Calculus
    import Optim
    using Optim: OptimizationTrace, OptimizationState, MultivariateOptimizationResults
    import NaNMath
    using Levenshtein

    export

        # distributions_ext.jl
        BetaAlt, GammaAlt,

        # settings.jl
        Setting, get_setting, default_settings!, default_test_settings!,

        # abstractdsgemodel.jl
        AbstractModel, transform_to_model_space!, transform_to_real_line!,
        n_states, n_shocks_exogenous, n_shocks_expectational, n_parameters, n_observables,
        spec, subspec, n_anticipated_shocks, n_states_augmented,
        n_presample_periods, zlb_start_index,
        dataroot, saveroot, inpath, workpath, rawpath, tablespath, figurespath, logpath,
        reoptimize, calculate_hessian,
        n_mh_blocks, n_mh_simulations, n_mh_burn, mh_thin, specify_mh_start, n_draws,
        data_vintage, cond_vintage, cond_id,
        specify_mode!, specify_hessian, load_parameters_from_file,
        use_population_forecast,


        # parameters.jl
        parameter, Transform, NullablePrior, AbstractParameter,
        Parameter, ParameterVector, ScaledParameter,
        UnscaledParameter, SteadyStateParameter, transform_to_real_line, transform_to_model_space,
        update, update!, transform_to_model_space, transform_to_real_line, Interval, ParamBoundsError,

        # statespace/
        Measurement, Transition, System, 

        # estimate/
        kalman_filter, kalman_filter_2part, likelihood, posterior, posterior!,
        optimize!, csminwel, hessian!, estimate, proposal_distribution,
        metropolis_hastings, compute_parameter_covariance, compute_moments,
        find_density_bands, prior,

        # forecast/
        kalman_smoother, disturbance_smoother, filter, compute_system,

        # models/
        steadystate!, Model990, SmetsWouters, eqcond, measurement,

        # solve/
        gensys, solve,

        # data/
        load_fred_data, load_data, transform_data, save_data_matrix!, hpfilter, difflog

    const VERBOSITY = Dict{Symbol,Int}(:none => 0, :low => 1, :high => 2)
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

end
