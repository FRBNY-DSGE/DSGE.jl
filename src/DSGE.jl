isdefined(Base, :__precompile__) && __precompile__()

module DSGE
    using Compat, Distributions, Roots.fzero, HDF5
    using DataStructures: SortedDict, insert!, ForwardOrdering
    using QuantEcon: solve_discrete_lyapunov
    
    if VERSION < v"0.4-"
        using Docile
    end

    export

        # distributions_ext.jl
        BetaAlt, GammaAlt,

        # settings.jl
        Setting, get_setting, default_settings, default_test_settings,

        # abstractdsgemodel.jl
        AbstractModel, transform_to_model_space!, transform_to_real_line!, n_states, n_shocks_exogenous,
        n_shocks_expectational, n_parameters, spec, subspec, 
        saveroot, inpath, workpath, rawpath, tablespath, figurespath, logpath,
        optimize, calculate_hessian, n_mh_blocks, n_mh_simulations,
        n_mh_burn, mh_thinning_step, data_vintage,

        
        # parameters.jl
        parameter, Transform, NullablePrior, AbstractParameter,
        Parameter, ParameterVector, ScaledParameter,
        UnscaledParameter, SteadyStateParameter, transform_to_real_line, transform_to_model_space,
        update, update!, transform_to_model_space, transform_to_real_line, Interval, ParamBoundsError,
        
        # estimate/
        kalman_filter, likelihood, posterior, posterior!,
        optimize!, csminwel, hessian!, estimate, proposal_distribution,
        metropolis_hastings, compute_parameter_covariance, compute_moments,
        make_moment_tables, find_density_bands, prior,
        
        # models/
        steadystate!, Model990, eqcond, measurement,

        # solve/
        ordschur, gensys, solve
    
    const VERBOSITY = @compat(Dict{Symbol,Int}(:none => 0, :low => 1, :high => 2))

    include("parameters.jl")
    include("distributions_ext.jl")
    include("abstractdsgemodel.jl")
    include("settings.jl")
    
    if VERSION <= v"0.4-"
        include("solve/ordered_qz.jl")
    end
    
    include("solve/gensys.jl")
    include("solve/solve.jl")
    
    include("estimate/kalman.jl")
    include("estimate/posterior.jl")
    include("estimate/optimize.jl")
    include("estimate/csminwel.jl")
    include("estimate/hessian.jl")
    include("estimate/estimate.jl")
    include("estimate/moments.jl")
    
    include("models/m990/m990.jl")
    include("models/m990/subspecs.jl")
    include("models/m990/eqcond.jl")
    include("models/m990/measurement.jl")
    
end
