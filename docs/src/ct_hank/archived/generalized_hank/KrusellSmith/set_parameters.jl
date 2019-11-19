function set_parameters()
    params = Dict{Symbol, Float64}()
    grids = Dict{Symbol, Any}()

    # Preferences
    params[:γ] = 2.        #coefficient of relative risk aversion
    params[:ρ] = 0.01         #rate of time preference

    # Production function
    params[:δ] = .025     # capital depreciation
    params[:α] = 1 / 3    # capital share

    # Aggregate shock
    params[:σ_tfp] = .007  # standard deviation of TFP shock
    params[:ρ_tfp] = .95     # quarterly autocorrelation of TFP shock

    # Idiosyncratic shocks
    params[:zz1] = 0.           # unemployed
    params[:zz2] = 1.           # employed
    grids[:z] = [params[:zz1],params[:zz2]]
    grids[:dz] = 1.

    # Transition probabilities
    params[:λ1] = 1 / 2  # expected duration of unemployment is 2 quarters
    params[:target_ur] = .07
    params[:target_er] = .93 # target employment rate
    params[:λ2] = (params[:λ1] / (params[:zz2] * params[:target_er] - params[:zz1]))*(params[:zz2] - params[:zz2] * params[:target_er]) # unemployment rate 7#
    grids[:lla] = [params[:λ1],params[:λ2]]

    # Tax system
    params[:μ] = .15        # UI replacement rate 15#
    params[:τ] = (params[:μ] / params[:zz2]) * (grids[:lla][2] /grids[:lla][1])	     # labor income tax

    # Labor supply
    grids[:zAvg] = (grids[:lla][1] * grids[:z][2] + grids[:lla][2] * grids[:z][1]) / (grids[:lla][1] + grids[:lla][2])

    ## Approximation Parameters
    # Wealth grid
    grids[:I] = 100
    grids[:amin] = 0.
    grids[:amax] = 100.
    grids[:a] = collect(linspace(grids[:amin],grids[:amax],grids[:I]))
    grids[:da] = (grids[:amax] - grids[:amin]) / (grids[:I] - 1)
    grids[:aa] = [grids[:a] grids[:a]]
    grids[:aaa] = reshape(grids[:aa], 2*grids[:I], 1)
    grids[:dz] = 1.
    grids[:J] = 2

    # Labor productivity
    grids[:zz] = ones(grids[:I],1) * grids[:z]'
    grids[:zzz] = reshape(grids[:zz], 2*grids[:I], 1)

    # Idiosyncratic shocks for income
    grids[:Aswitch] = [-speye(grids[:I]) * grids[:lla][1] speye(grids[:I]) * grids[:lla][1]; speye(grids[:I]) * grids[:lla][2] -speye(grids[:I]) * grids[:lla][2]]

    init_params = Dict{Symbol, Float64}()
    approx_params = Dict{Symbol, Any}()

    # Steady state computations
    init_params[:rmin] = .0001       # lower bound for steady state interest rate
    init_params[:rmax] = params[:ρ]        # upper bound for steady state interest rate
    init_params[:r0] = .005          # initial guess for steady state interest rate
    approx_params[:maxit_HJB] = 100        # maximum iterations on steady state HJB
    approx_params[:crit] = 1e-6
    approx_params[:Δ] = 1e4
    approx_params[:crit_HJB] = 1e-6        # error criterion for steady state value function convergence
    approx_params[:Δ_HJB] = 1e4        # update size for implicit scheme on steady state HJ
    approx_params[:Ir] = 100           # maximum iterations on steady state interest rate
    approx_params[:crit_S] = 1e-5     # error criterion for steady state interest rate

    # Number of variables in the system
    grids[:n_jump_vars] = 2 * grids[:I]
    grids[:n_state_vars] = 2 * grids[:I]-1 + 1
    grids[:n_static_conditions] = 6
    grids[:nVars] = grids[:n_jump_vars] + grids[:n_state_vars] + grids[:n_static_conditions]
    grids[:nEErrors] = 2 * grids[:I]

    # Krylov and value function reduction parameters
    n_knots = 12
    c_power = 7
    knots = collect(linspace(grids[:amin], grids[:amax], n_knots - 1))
    knots = (grids[:amax] - grids[:amin])/(2^c_power - 1) * ((knots - grids[:amin]) / (grids[:amax] - grids[:amin]) + 1).^c_power + grids[:amin] - (grids[:amax] - grids[:amin])/(2^c_power - 1)
    krylov_dim = 5

    # Create ReductionData variable: reduceDistribution, reduceV, krylov_dim,
    red_params = ReductionData(true, true, krylov_dim, Dict(1 => knots), grids[:a], 1, 2)

    return params, init_params, grids, approx_params, red_params
end

