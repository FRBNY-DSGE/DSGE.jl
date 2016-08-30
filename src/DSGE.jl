isdefined(Base, :__precompile__) && __precompile__()

module DSGE
    using Distributions, HDF5
    using DataStructures: SortedDict, insert!, ForwardOrdering, OrderedDict
    using FredData, DataFrames, Base.Dates
    using QuantEcon: solve_discrete_lyapunov
    import Calculus
    using Roots: fzero, ConvergenceFailed
    import Optim: optimize, Optimizer
    # using Optim: OptimizationTrace, OptimizationState, MultivariateOptimizationResults
    # using Optim: OptimizationOptions, optimize, SimulatedAnnealing, Optimizer
    using Optim

    export

        # distributions_ext.jl
        BetaAlt, GammaAlt, RootInverseGamma,

        # settings.jl
        Setting, get_setting, default_settings!, default_test_settings!,

        # abstractdsgemodel.jl
        AbstractModel, transform_to_model_space!, transform_to_real_line!,
        description,
        n_states, n_states_augmented, n_shocks_exogenous, n_shocks_expectational,
        n_equilibrium_conditions, n_observables, n_parameters, n_parameters_steady_state,
        n_parameters_free, n_anticipated_shocks,
        spec, subspec,
        dataroot, saveroot, inpath, workpath, rawpath, tablespath, figurespath, logpath,
        reoptimize, calculate_hessian,
        n_mh_blocks, n_mh_simulations, n_mh_burn, mh_thin, specify_mh_start,
        data_vintage,
        specify_mode!, specify_hessian, load_parameters_from_file,
        use_population_forecast,

        # parameters.jl
        parameter, Transform, NullablePrior, AbstractParameter,
        Parameter, ParameterVector, ScaledParameter,
        UnscaledParameter, SteadyStateParameter, transform_to_real_line, transform_to_model_space,
        update, update!, transform_to_model_space, transform_to_real_line, Interval, ParamBoundsError,

        # estimate/
        kalman_filter, likelihood, posterior, posterior!,
        optimize!, csminwel, simulated_annealing, hessian!, estimate, proposal_distribution,
        metropolis_hastings, compute_parameter_covariance, compute_moments,
        find_density_bands, prior,

        # models/
        steadystate!, Model990, SmetsWouters, eqcond, measurement,

        # solve/
        gensys, solve,

        # data/
        load_data, load_data_levels, load_fred_data, transform_data, save_data,
        df_to_matrix, hpfilter, difflog, quartertodate, percapita, nominal_to_real,
        hpadjust, oneqtrpctchange, annualtoquarter

    const VERBOSITY = Dict(:none => 0, :low => 1, :high => 2)
    const DSGE_DATE_FORMAT = "yymmdd"

    include("parameters.jl")
    include("distributions_ext.jl")
    include("abstractdsgemodel.jl")
    include("settings.jl")
    include("defaults.jl")
    include("util.jl")

    include("data/load_data.jl")
    include("data/fred_data.jl")
    include("data/transformations.jl")
    include("data/transform_data.jl")
    include("data/util.jl")

    include("solve/gensys.jl")
    include("solve/solve.jl")

    include("estimate/kalman.jl")
    include("estimate/posterior.jl")
    include("estimate/optimize.jl")
    include("estimate/csminwel.jl")
    include("estimate/hessian.jl")
    include("estimate/hessizero.jl")
    include("estimate/simulated_annealing.jl")
    include("estimate/estimate.jl")
    include("estimate/moments.jl")

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
