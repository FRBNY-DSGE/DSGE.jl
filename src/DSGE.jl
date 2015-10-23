module DSGE
    using Compat
    using Distributions, Roots.fzero, MATLAB, HDF5

    export
        # DSGE.jl
        savepath, inpath, outpath, tablepath, plotpath, logpath,

        # distributions_ext.jl
        PointMass, BetaAlt, GammaAlt,

        # abstractdsgemodel.jl
        AbstractDSGEModel, ModelInds, makedict, tomodel!, toreal!, #toreal, tomodel, update!,  Parameters,

        # new_parameters.jl
        parameter, Transform, Untransformed, SquareRoot, Exponential, NullablePrior,
        AbstractParameter, Parameter, ParameterVector, ScaledParameter, UnscaledParameter, SteadyStateParameter,
        toreal, tomodel, update, update!, tomodel, toreal, Interval,

        # solve/
        ordschur, gensys, solve,

        # estimate/
        dlyap!, kalcvf2NaN, kalsmth_k93, prior, likelihood, posterior, posterior!, csminwel, hessizero!, estimate, proposal_distribution, metropolis_hastings,

        # models/
        steadystate!, Model990, model_specifications, eqcond, measurement

    include("new_parameters.jl")
    include("distributions_ext.jl")
    include("abstractdsgemodel.jl")
    # include("parameters.jl")


    if VERSION <= v"0.4-"
        include("solve/ordered_qz.jl")
    end

    include("solve/gensys.jl")
    include("solve/solve.jl")

    include("estimate/kalman.jl")
    include("estimate/posterior.jl")
    include("estimate/csminwel.jl")
    include("estimate/hessian.jl")
    include("estimate/estimate.jl")

    include("models/m990/m990.jl")
    include("models/m990/eqcond.jl")
    include("models/m990/measurement.jl")
end
