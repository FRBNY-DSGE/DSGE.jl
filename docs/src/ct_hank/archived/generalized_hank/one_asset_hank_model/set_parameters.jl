function set_parameters(n_v, n_g, n_p, n_shocks; n_exp_errors::Int64 = n_v)
    #######################
    # Economic Parameters
    #######################

    # utility function
    params = Dict{Symbol, Float64}()

    params[:coefrra]    = 1.0
    params[:frisch]     = 0.5
    params[:meanlabeff] = 3.0
    params[:maxhours]   = 1.0

    params[:ceselast]          = 10.
    params[:priceadjust]       = 100.
    params[:taylor_inflation]  = 1.25  # taylor rule coefficient on inflation
    params[:taylor_outputgap]  = 0.    # taylor rule coefficient on output
    params[:labtax]            = 0.2   # marginal tax rate
    params[:govbondtarget]     = 6.	   # multiple of quarterly GDP
    params[:lumptransferpc]    = 0.06  # 6% of quarterly GDP in steady state
    params[:govbcrule_fixnomB] = 0.

    params[:labdisutil] = params[:meanlabeff] / ( (0.75 ^(-params[:coefrra])) * ((1./3.)^(1/params[:frisch])))

    # Aggregate shocks
    params[:ssigma_MP] = sqrt(0.05)
    params[:ttheta_MP] = 0.25

    #####################################
    # Initial Guess Parameters
    #####################################

    init_params = Dict{Symbol, Float64}()

    # steady state r
    init_params[:r0]   = 0.005
    init_params[:rmax] = 0.08
    init_params[:rmin] = 0.001

    # steady state rho iteration
    init_params[:rrho0]  = 0.02
    init_params[:rhomax] = 0.05
    init_params[:rhomin] = 0.005

    ####################################
    # Set Grids and relevant parameters
    ####################################

    grids = Dict{Symbol, Any}()

    #################
    # Asset
    #################
    I          = 100 # Number of grid points
    agridparam = 1   # Bending coefficient - 1 for linear
    amin       = 0.  # Minimum asset grid value
    amax       = 40. # Maximum asset grid value

    grids[:I] = I
    params[:amin] = amin
    params[:amax] = amax
    grids[:a] = construct_asset_grid(I, agridparam, amin, amax)

    #################
    # Income
    #################
    J = 2                                           # Number of income states
    ygrid_combined = [0.2, 1]                       # Raw levels of income earned in each income state
    ymarkov_combined = [-0.5 0.5; 0.0376 -0.0376]   # Markov transition matrix from income states

    g_z = compute_stationary_income_distribution(ymarkov_combined, J)
    zz  = construct_labor_income_grid(ygrid_combined, g_z, params[:meanlabeff], I)

    grids[:J] = J
    grids[:g_z] = g_z
    grids[:zz] = zz
    grids[:z] = zz[1,:]
    grids[:ymarkov_combined] = ymarkov_combined

    ###################
    # Grid Dimensions
    ###################
    grids[:n_jump_vars] = n_v         # jump variables and inflation
    grids[:n_state_vars] = n_g             # endogenous state variables
    grids[:n_static_conditions] = n_p                    # Number of static relations via market clearing
    grids[:n_shocks] = n_shocks                # Monetary policy only
    grids[:n_exp_errors] = n_exp_errors
    grids[:n_vars] = n_v + n_g + n_p

    ##############################################
    # Approximation Parameters for Steady State
    ##############################################
    approx_params = Dict{Symbol, Any}()

    approx_params[:Ir] = 100 # number of iterations to find v, g, and p fixed points
    approx_params[:maxit_HJB] = 500
    approx_params[:tol_HJB] = 1e-8
    approx_params[:Delta_HJB] = 1e6
    approx_params[:maxit_KFE] = 1000
    approx_params[:tol_KFE] = 1e-12
    approx_params[:Delta_KFE] = 1e6
    approx_params[:niter_hours] = 10
    approx_params[:IterateR] = false
    approx_params[:IterateRho] = true
    approx_params[:crit_S] = 1e-5

    #######################################
    # Approximation Grid for ValueFunction
    #######################################
    n_knots = 12
    c_power = 1
    ss_array = grids[:a]
    n_post = length(zz[1,:])
    n_prior = 1
    general_spline = true
    knots = collect(linspace(amin, amax, n_knots - 1))
    knots1 = (amax-amin)/(2^c_power-1)*((knots-amin)/(amax-amin)+1).^c_power + amin - (amax-amin)/(2^c_power-1)
    knots_dict = Dict(1 => knots1)
    krylov_dim = 20
    reduce_distribution = true
    reduce_v = true
    red_params = ReductionData(reduce_distribution, reduce_v, krylov_dim, knots_dict, ss_array, n_prior, n_post)

    return params, init_params, grids, approx_params, red_params
end