function set_parameters()

    # This function sets parameters

    ######################
    # PARAMETERS
    ######################
    params = Dict{Symbol, Float64}()

    params[:ga] = 2.       # CRRA utility with parameter gamma
    params[:alpha] = 0.35; # Production function F = K^alpha * L^(1-alpha)
    params[:delta] = 0.1;  # Capital depreciation
    params[:zmean] = 1.0;      # mean O-U process (in levels). This parameter has to be adjusted to ensure that the mean of z (truncated gaussian) is 1.
    params[:sig2] = (0.10)^2;  # sigma^2 O-U
    params[:Corr] = exp(-0.3);  # persistence -log(Corr)  O-U
    params[:rho] = 0.05;   # discount rate
    params[:K] = 3.8;      # initial aggregate capital. It is important to guess a value close to the solution for the algorithm to converge

    ######################
    # Initial Guesses
    ######################
    init_params = Dict{Symbol, Float64}

    init_params[:r] = params[:alpha] * params[:K]^(params[:alpha] - 1) - params[:delta] # interest rate
    init_params[:w] = (1 - params[:alpha]) * params[:K]^params[:alpha]                  # wage

    ####################
    # Grids
    ####################
    grids = Dict{Symbol, Any}()

    grids[:J] = 40;         # number of z points
    grids[:zmin] = 0.5;   # Range z
    grids[:zmax] = 1.5;
    grids[:amin] = -1.    # borrowing constraint
    grids[:amax] = 30.;    # range a
    grids[:I] = 100;        # number of a points
    grids[:a] = collect(linspace(amin, amax, I)) # wealth vector
    grids[:da] = (grids[:amax] - grids[:amin])/(grids[:I] - 1)
    grids[:z] = collect(linspace(grids[:zmin], grids[:zmax], grids[:J]))'
    grids[:dz] = (grids[:zmax] - grids[:zmin])/(grids[:J] - 1)
    grids[:dz2] = dz*dz
    grids[:aa] = grids[:a]*ones(1, grids[:J])
    grids[:zz] = ones(grids[:I], 1)*grids[:z]

    # ORNSTEIN-UHLENBECK IN LEVELS: defining the income process
    the = params[:-log(Corr)]
    Var = params[:sig2]/(2*the)
    grids[:mu] = the*(params[:zmean] - grids[:z])
    grids[:s2] = params[:sig2].*ones(1, grids[:J])

    ##############################
    # Approximation Parameters
    ##############################
    approx_params = Dict{Symbol, Any}()

    approx_params[:relax] = 0.99; # relaxation parameter
    approx_params[:maxit_HJB]  = 100.     # maximum number of iterations in the HJB loop
    approx_params[:maxit_KFE] = 100;    # maximum number of iterations in the KFE loop
    approx_params[:HJB_tol] = 10^(-6.); # tolerance HJB loop
    approx_params[:KFE_tol] = 10^(-5.);   # tolerance KFE loop
    approx_params[:HJB_Delta] = 1000.;   # Delta in HJB algorithm

    ##############################
    # Value Function Reduction
    ##############################
    approx_valuef = Dict{Symbol, Any}

    return params, init_params, grids, approx_params, approx_valuef
end