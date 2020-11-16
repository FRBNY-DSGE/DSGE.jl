"""
```
TwoAssetHANK{T} <: AbstractCTModel{T}
```

The `TwoAssetHANK` type defines the structure of the two-asset HANK model
originally written in MATLAB by Ben Moll and Se Hyoun Ahn.

### Fields

#### Parameters and Steady-States
* `parameters::Vector{AbstractParameter}`: Vector of all time-invariant model
  parameters.

* `steady_state::Vector{AbstractParameter}`: Model steady-state values, computed
  as a function of elements of `parameters`.

* `keys::OrderedDict{Symbol,Int}`: Maps human-readable names for all model
  parameters and steady-states to their indices in `parameters` and
  `steady_state`.

#### Inputs to Measurement and Equilibrium Condition Equations

The following fields are dictionaries that map human-readable names to row and
column indices in the matrix representations of of the measurement equation and
equilibrium conditions.

* `endogenous_states::OrderedDict{Symbol,Int}`: Maps each state to a column in
  the measurement and equilibrium condition matrices.

* `exogenous_shocks::OrderedDict{Symbol,Int}`: Maps each shock to a column in
  the measurement and equilibrium condition matrices.

* `expected_shocks::OrderedDict{Symbol,Int}`: Maps each expected shock to a
  column in the measurement and equilibrium condition matrices.

* `equilibrium_conditions::OrderedDict{Symbol,Int}`: Maps each equlibrium
  condition to a row in the model's equilibrium condition matrices.

* `endogenous_states_augmented::OrderedDict{Symbol,Int}`: Maps lagged states to
  their columns in the measurement and equilibrium condition equations. These
  are added after `gensys` solves the model.

* `observables::OrderedDict{Symbol,Int}`: Maps each observable to a row in the
  model's measurement equation matrices.

* `pseudo_observables::OrderedDict{Symbol,Int}`: Maps each pseudo-observable to
  a row in the model's pseudo-measurement equation matrices.

#### Model Specifications and Settings

* `spec::String`: The model specification identifier, \"one_asset_hank\", cached
  here for filepath computation.

* `subspec::String`: The model subspecification number, indicating that some
  parameters from the original model spec (\"ss0\") are initialized
  differently. Cached here for filepath computation.

* `settings::Dict{Symbol,Setting}`: Settings/flags that affect computation
  without changing the economic or mathematical setup of the model.

* `test_settings::Dict{Symbol,Setting}`: Settings/flags for testing mode

#### Other Fields

* `rng::MersenneTwister`: Random number generator. Can be is seeded to ensure
  reproducibility in algorithms that involve randomness (such as
  Metropolis-Hastings).

* `testing::Bool`: Indicates whether the model is in testing mode. If `true`,
  settings from `m.test_settings` are used in place of those in `m.settings`.

* `observable_mappings::OrderedDict{Symbol,Observable}`: A dictionary that
  stores data sources, series mnemonics, and transformations to/from model units.
  DSGE.jl will fetch data from the Federal Reserve Bank of
  St. Louis's FRED database; all other data must be downloaded by the
  user. See `load_data` and `Observable` for further details.

* `pseudo_observable_mappings::OrderedDict{Symbol,PseudoObservable}`: A
  dictionary that stores names and transformations to/from model units. See
  `PseudoObservable` for further details.
"""
mutable struct TwoAssetHANK{T} <: AbstractCTModel{T}
    parameters::ParameterVector{T}                          # vector of all time-invariant model parameters
    steady_state::ParameterVector{T}                        # model steady-state values
    keys::OrderedDict{Symbol,Int}                           # human-readable names for all the model
                                                            # parameters and steady-states

    endogenous_states::OrderedDict{Symbol,Vector{Int}}      # these fields used to create matrices in the
    exogenous_shocks::OrderedDict{Symbol,Int}               # measurement & equilibrium condition equations
    expected_shocks::OrderedDict{Symbol,Vector{Int}}        #
    equilibrium_conditions::OrderedDict{Symbol,Vector{Int}} #
    endogenous_states_augmented::OrderedDict{Symbol,Int}    #
    observables::OrderedDict{Symbol,Int}                    #
    pseudo_observables::OrderedDict{Symbol,Int}             #

    spec::String                                            # Model specification number (eg "m990")
    subspec::String                                         # Model subspecification (eg "ss0")
    settings::Dict{Symbol,Setting}                          # Settings/flags for computation
    test_settings::Dict{Symbol,Setting}                     # Settings/flags for testing mode
    rng::MersenneTwister                                    # Random number generator
    testing::Bool                                           # Whether we are in testing mode or not

    observable_mappings::OrderedDict{Symbol, Observable}
    pseudo_observable_mappings::OrderedDict{Symbol, PseudoObservable}
end

description(m::TwoAssetHANK) = "Julia implementation of the two-asset HANK model
originally written in MATLAB by Ben Moll and Se Hyoun Ahn."


"""
`init_model_indices!(m::TwoAssetHANK)`

Arguments:
`m:: TwoAssetHANK`: a model object

Description:
Initializes indices for all of `m`'s states, shocks, and equilibrium conditions.
"""
function init_model_indices!(m::TwoAssetHANK)

    n_v = get_setting(m, :n_v)
    n_g = get_setting(m, :n_g)
    n_p = get_setting(m, :n_p)
    n_Z = get_setting(m, :n_Z)

    # Endogenous states
    endogenous_states = collect([:value_function, :distribution, :K, :r_b, :Y, :C, :Z])

    # Exogenous shocks
    exogenous_shocks = collect([:mp_shock])

    # Expectations shocks
    expected_shocks = collect([:E_V, :E_π])

    # Equilibrium conditions
    equilibrium_conditions = collect([:value_function, :distribution, :K, :r_b, :Y, :C, :Z])

    # Additional states added after solving model
    # Lagged states and observables measurement error
    endogenous_states_augmented = []

    # Observables
    observables = keys(m.observable_mappings)

    # Pseudo-observables
    pseudo_observables = keys(m.pseudo_observable_mappings)

    # Assign dimensions
    m.endogenous_states[:value_function]  = collect(1:n_v)
    m.endogenous_states[:distribution]    = collect(n_v + 1:n_v + n_g)
    m.endogenous_states[:K]               = [n_v + n_g + 1]
    m.endogenous_states[:r_b]             = [n_v + n_g + 2]

    m.endogenous_states[:Y]               = [n_v + n_g + 3]
    m.endogenous_states[:C]               = [n_v + n_g + 4]

    m.endogenous_states[:Z]               = [n_v + n_g + n_p + 1]

    m.expected_shocks[:E_V]               = m.endogenous_states[:value_function]
    #m.expected_shocks[:E_π]               = m.endogenous_states[:inflation]

    for (i,k) in enumerate(endogenous_states_augmented)
        m.endogenous_states_augmented[k] = i + length(endogenous_states)
    end
    for (i,k) in enumerate(observables);          m.observables[k]          = i end
    for (i,k) in enumerate(exogenous_shocks);     m.exogenous_shocks[k]     = i end
    for (i,k) in enumerate(pseudo_observables);   m.pseudo_observables[k]   = i end
end


function TwoAssetHANK(subspec::String="ss0";
                      custom_settings::Dict{Symbol, Setting} = Dict{Symbol, Setting}(),
                      testing = false)

    # Model-specific specifications
    spec               = split(basename(@__FILE__),'.')[1]
    subspec            = subspec
    settings           = Dict{Symbol,Setting}()
    test_settings      = Dict{Symbol,Setting}()
    rng                = MersenneTwister(0)

    # Initialize empty model
    m = TwoAssetHANK{Float64}(
            # Model parameters and steady state values # before 2nd arg was Vector{Float64}()
            Vector{AbstractParameter{Float64}}(), Vector{Float64}(),
            OrderedDict{Symbol,Int}(),

            # Model indices
            OrderedDict{Symbol,Vector{Int}}(), OrderedDict{Symbol,Int}(),
            OrderedDict{Symbol,Vector{Int}}(), OrderedDict{Symbol,Vector{Int}}(),
            OrderedDict{Symbol,Int}(),         OrderedDict{Symbol,Int}(),
            OrderedDict{Symbol,Int}(),

            spec,
            subspec,
            settings,
            test_settings,
            rng,
            testing,
            OrderedDict{Symbol,Observable}(),
            OrderedDict{Symbol,PseudoObservable}())

    # Initialize parameters
    init_parameters!(m)

    # Set settings
    model_settings!(m)
    default_test_settings!(m)
    adjust_parameters!(m)

    init_model_indices!(m)

    for custom_setting in values(custom_settings)
        m <= custom_setting
    end

    # Test to see if necessary values are defined, throw error otherwise
#=    try
        @assert typeof(get_setting(m, :reduce_state_vars)) == Bool
        @assert typeof(get_setting(m, :reduce_v)) == Bool
    catch
        error("Need to specify in settings whether to perform distribution
               and/or value function reduction")
    end

    if get_setting(m, :reduce_state_vars) || get_setting(m, :reduce_v)
        try
            @assert get_setting(m, :n_jump_vars) > 0
            @assert get_setting(m, :n_state_vars) > 0
            @assert get_setting(m, :n_state_vars_unreduce) >= 0
        catch
            error("Need to enter values in settings for n_jump_vars, n_state_vars,
                   and/or n_state_vars_unreduce in Settings. These indicate the
                   number of jump and aggregate state variables.")
        end
    end

    if get_setting(m, :reduce_state_vars)
        try
            @assert get_setting(m, :krylov_dim) > 0
        catch
            error("Need to specify in settings desired dimension of
                   Krylov subspace, krylov_dim.")
        end
    end

    if get_setting(m, :reduce_v)
        try
            @assert get_setting(m, :n_prior) > 0 && typeof(get_setting(m, :n_prior)) == Int
            @assert get_setting(m, :n_post) > 0 && typeof(get_setting(m, :n_post)) == Int
            @assert in(:knots_dict, keys(m.settings))
            @assert in(:spline_grid, keys(m.settings))
        catch
            error("Need to specify in settings n_prior, n_post, knots_dict
                   and/or spline_grid for value function reduction.
                   See reduction.jl for their uses.")
        end
    end
=#
    # Compute steady state
#    steadystate!(m)
    return m
end

"""
```
init_parameters!(m::TwoAssetHANK)
```

Initializes the model's parameters, as well as empty values for the steady-state
parameters (in preparation for `steadystate!(m)` being called to initialize
those).
"""
function init_parameters!(m::TwoAssetHANK)
    # Preferences
    m <= parameter(:coefrra, 1.0, fixed=true,
                   description="Relative risk aversion",
                   tex_label="coefrra")

    m <= parameter(:rrho, 0.014526858925819, fixed=true,
                   description="Discount factor (quarterly)",
                   tex_label="rrho")

    m <= parameter(:ddeath, 1.0 / (4.0 * 45.0), fixed=true,
                   description="Death rate (quarterly)",
                   tex_label="ddeath")

    m <= parameter(:ggamma, 1.0, fixed=true,
                   description="Risk aversion",
                   tex_label="ggama")

    # Technology
    m <= parameter(:aalpha, 0.4, fixed=true,
                   description="Capital share",
                   tex_label="aalpha")

    m <= parameter(:nnu_aggZ, 0.05, fixed=true,
                   description="nnu_aggZ",
                   tex_label="nnu_aggZ")

    m <= parameter(:ssigma_aggZ, 0.007 / (1.0 - m[:aalpha]), fixed=true,
                   description="ssigma_aggZ",
                   tex_label="ssigma_aggZ")

    m <= parameter(:ddelta, 0.075 / 4.0 , fixed=true,
                   description="Depreciation rate (quarterly)",
                   tex_label="ddelta")

    # Portfolio adjustment cost function
    m <= parameter(:chi0, 0.008198815318426, fixed=true,
                   description="Portfolio adjustment cost function",
                   tex_label="chi0")

    m <= parameter(:chi1, 0.072974260992256, fixed=true,
                   description="Portfolio adjustment cost function",
                   tex_label="chi1")

    m <= parameter(:chi2, 1.3583017016896, fixed=true,
                   description="Portfolio adjustment cost function",
                   tex_label="chi2")

    m <= parameter(:a_lb, 0.02, fixed=true,
                   description="Portfolio adjustment cost function",
                   tex_label="a_lb")

    # Tax code
    m <= parameter(:trans, 0.1, fixed=true,
                   description="Transfer to households",
                   tex_label="trans")

    m <= parameter(:tau_I, 0.30, fixed=true,
                   description="Tax rate on wage income",
                   tex_label="tau_I")

    # Deposit share
    m <= parameter(:xxi, 0.0, fixed=true,
                   description="Deposit share",
                   tex_label="xxi")

    # Perfect annuity markets
    m <= parameter(:pam, 1.0, fixed=true,
                   description="Perfect annuity markets",
                   tex_label="pam")

    # Liquid return
    m <= parameter(:r_b_SS, 0.005, fixed=true,
                   description="Steady-state savings rate, liquid assets",
                   tex_label="r_b_SS")
    m <= parameter(:borrwedge_SS, 0.020137162949122, fixed=true,
                   description="Steady-state borrowing wedge, liquid assets",
                   tex_label="borrwedge_SS")
    m <= parameter(:r_b_borr_SS, m[:r_b_SS] + m[:borrwedge_SS], fixed=true,
                   description="Steady-state borrowing rate, liquid assets",
                   tex_label="r_b_borr_SS")
    m <= parameter(:pphi, 2.0, fixed=true,
                   description="Supply elasticity (only relevant when r_b_phi = 1)",
                   tex_label="pphi")

    # Stochastic Properties of aggregate TFP shock
    m <= parameter(:nnu_aggZ, 0.25, fixed=true,
                   description="Persistence parameter in Ornstein-Uhlenbeck process",
                   tex_label="nnu_aggZ")
    m <= parameter(:sigma_aggZ, 0.007 / (1 - m[:aalpha]), fixed=true,
                   description="Volatility parameter in Ornstein-Uhlenbeck process",
                   tex_label="sigma_aggZ")

    # Household effects of TFP shock
    m <= parameter(:kappa, 10.0, fixed=true,
                   description="Household effects of TFP shock",
                   tex_label="kappa")

    # Steady State Parameters
    m <= SteadyStateParameter(:KL_SS, NaN,
                              description = "")
    m <= SteadyStateParameter(:r_a_SS, NaN,
                              description = "")
    m <= SteadyStateParameter(:w_SS, NaN,
                              description = "")
    m <= SteadyStateParameter(:K_SS, NaN,
                              description = "")
    m <= SteadyStateParameterGrid(:u_SS, Vector{Float64}(),
                                  description = "")
    m <= SteadyStateParameterGrid(:c_SS, Vector{Float64}(),
                                  description = "")
    m <= SteadyStateParameterGrid(:d_SS, Vector{Float64}(),
                                  description = "")
    m <= SteadyStateParameterGrid(:V_SS, Vector{Float64}(),
                                  description = "")
    m <= SteadyStateParameterGrid(:g_SS, Vector{Float64}(),
                                  description = "")
    m <= SteadyStateParameterGrid(:dab, Vector{Float64}(),
                                  description = "")
    m <= SteadyStateParameterGrid(:dab_g, Vector{Float64}(),
                                  description = "")
m <= SteadyStateParameter(:C_SS, NaN,
                          description = "")
m <= SteadyStateParameter(:C_Var_SS, NaN,
                          description = "")
m <= SteadyStateParameter(:I_SS, NaN,
                          description = "")
m <= SteadyStateParameter(:B_SS, NaN,
                          description = "")
m <= SteadyStateParameter(:Y_SS, NaN,
                          description = "")
m <= SteadyStateParameter(:n_SS, NaN,
                          description = "")
m <= SteadyStateParameterGrid(:earn_SS, Vector{Float64}(),
                              description = "")
m <= SteadyStateParameter(:earn_Var_SS, NaN,
                          description = "")
m <= SteadyStateParameterGrid(:IcF_SS, Vector{Float64}(),
                              description = "")
m <= SteadyStateParameterGrid(:IcB_SS, Vector{Float64}(),
                              description = "")
m <= SteadyStateParameterGrid(:Ic0_SS, Vector{Float64}(),
                              description = "")
m <= SteadyStateParameterGrid(:IcFB_SS, Vector{Float64}(),
                              description = "")
m <= SteadyStateParameterGrid(:IcBF_SS, Vector{Float64}(),
                              description = "")
m <= SteadyStateParameterGrid(:IcBB_SS, Vector{Float64}(),
                              description = "")
m <= SteadyStateParameterGrid(:Ic00_SS, Vector{Float64}(),
                              description = "")
m <= SteadyStateParameterGrid(:IcBB_SS, Vector{Float64}(),
                              description = "")
m <= SteadyStateParameterGrid(:WHTM_indicator, Vector{Float64}(),
                              description = "")
m <= SteadyStateParameter(:WHTM_SS, NaN,
                          description = "")
m <= SteadyStateParameter(:C_WHTM_SS, NaN,
                          description = "")
m <= SteadyStateParameterGrid(:PHTM_indicator, Vector{Float64}(),
                              description = "")
m <= SteadyStateParameter(:PHTM_SS, NaN,
                          description = "")
m <= SteadyStateParameter(:C_PHTM_SS, NaN,
                          description = "")
m <= SteadyStateParameterGrid(:NHTM_indicator, Vector{Float64}(),
                              description = "")
m <= SteadyStateParameter(:NHTM_SS, NaN,
                          description = "")
m <= SteadyStateParameter(:C_NHTM_SS, NaN,
                          description = "")
m <= SteadyStateParameterGrid(:gg_SS, Vector{Float64}(),
                              description = "")
m <= SteadyStateParameter(:illiquid_wedge, NaN,
                          description = "")
m <= SteadyStateParameterGrid(:vars_SS, Vector{Float64}(),
                              description = "")
end

"""
```
steadystate!(m::TwoAssetHANK)
```

Calculates the model's steady-state values. `steadystate!(m)` must be called whenever
the parameters of `m` are updated.
"""
function steadystate!(m::TwoAssetHANK)
    # Read in parameters
    aalpha = m[:aalpha].value
    ddelta = m[:ddelta].value
    ddeath = m[:ddeath].value
    rrho   = m[:rrho].value
    chi0   = m[:chi0].value
    chi1   = m[:chi1].value
    chi2   = m[:chi2].value
    a_lb   = m[:a_lb].value
    pam    = m[:pam].value
    xxi    = m[:xxi].value
    ggamma = m[:ggamma].value
    tau_I  = m[:tau_I].value
    trans  = m[:trans].value

    # Set liquid rates
    r_b_SS       = m[:r_b_SS].value
    r_b_borr_SS  = m[:r_b_borr_SS].value
    borrwedge_SS = m[:borrwedge_SS].value

    lambda   = get_setting(m, :lambda)
    K_liquid = get_setting(m, :K_liquid)

    aggregate_variables        = get_setting(m, :aggregate_variables)
    distributional_variables   = get_setting(m, :distributional_variables)
    distributional_variables_1 = get_setting(m, :distributional_variables_1)

    # Read in approximation parameters
    maxit_HJB  = get_setting(m, :maxit_HJB)
    crit_HJB   = get_setting(m, :crit_HJB)
    Delta      = get_setting(m, :Delta)

    maxit_HIS  = get_setting(m, :maxit_HIS)
    crit_HIS   = get_setting(m, :crit_HIS)
    start_HIS  = get_setting(m, :start_HIS)

    maxit_KFE  = get_setting(m, :maxit_KFE)
    crit_KFE   = get_setting(m, :crit_KFE)
    Delta_KFE  = get_setting(m, :Delta_KFE)

    maxit_KL   = get_setting(m, :maxit_KL)
    crit_KL    = get_setting(m, :crit_KL)
    relax_KL   = get_setting(m, :relax_KL)

    # Read in grids
    I       = get_setting(m, :I)
    J       = get_setting(m, :J)
    a_g     = get_setting(m, :a_g)
    b_g     = get_setting(m, :b_g)
    I_g     = get_setting(m, :I_g)
    J_g     = get_setting(m, :J_g)
    N       = get_setting(m, :N)
    a       = get_setting(m, :a)
    b       = get_setting(m, :b)
    y       = get_setting(m, :y)
    y_dist  = get_setting(m, :y_dist)
    y_mean  = get_setting(m, :y_mean)
    KL      = get_setting(m, :KL_0)

    interp_decision = get_setting(m, :interp_decision)

    # Set liquid rates
    r_b      = r_b_SS
    r_b_borr = r_b_borr_SS

    # Compute prices associated with initial guess of KL
    w	= (1 - aalpha) * (KL ^ aalpha)
    r_a	= aalpha * (KL ^ (aalpha - 1)) - ddelta

    a_grid, a_g_grid, b_grid, b_g_grid, y_grid, y_g_grid, r_b_grid, r_b_g_grid, daf_grid, daf_g_grid, dab_grid, dab_g_grid, dab_tilde_grid, dab_g_tilde_grid, dab_g_tilde_mat, dab_g_tilde, dbf_grid, dbf_g_grid, dbb_grid, dbb_g_grid = set_grids(a, b, a_g, b_g, vec(y), r_b, r_b_borr)
    r_b_vec   = r_b .* (b   .>= 0) + r_b_borr .* (b   .< 0)

    # Construct problem functions
    util, deposit, cost = construct_problem_functions(ggamma, chi0, chi1, chi2, a_lb)

    # Initial consumption and value function
    c_0 = (1-xxi) * w * y_grid .+ (r_a + ddeath*pam) .* a_grid +
        (r_b_borr + ddeath * pam) .* b_grid .+ trans - tau_I * w * y_grid
    V_0	= 1/(rrho + ddeath) * util.(c_0)

    # Initial distribution
    gg0 = zeros(I_g, J_g, N)
    gg0[b .== 0, 1, :] = vec(y_dist)
    gg0 = gg0 ./ sum(gg0)
    gg0 = gg0 ./ dab_g_tilde_grid 	# ensure integration to 1
    gg0 = vec(gg0) #, I_g*J_g*N, 1)
    gg = gg0

    aau = Array{Float64}(undef, 0, 0)
    bbu = Array{Float64}(undef, 0, 0)
    ccu = Array{Float64}(undef, 0, 0)

    Vn1 = Array{Float64}(undef, I, J, N)
    Vn  = Array{Float64}(undef, I, J, N)

    gg_tilde = Array{Float64}(undef, 0, 0)

    g = Array{Float64}(undef, 0, 0)
    c = Array{Float64}(undef, I, J, N)
    s = Array{Float64}(undef, I, J, N)
    u = Array{Float64}(undef, I, J, N)
    d = Array{Float64}(undef, I, J, N)
    c_g = Array{Float64}(undef, 0, 0)
    s_g = Array{Float64}(undef, 0, 0)
    d_g = Array{Float64}(undef, 0, 0)

    IcF = BitArray{3}(undef, I,J,N)
    IcB = BitArray{3}(undef, I,J,N)
    Ic0 = BitArray{3}(undef, I,J,N)
    IcFB = BitArray{3}(undef, I,J,N)
    IcBF = BitArray{3}(undef, I,J,N)
    IcBB = BitArray{3}(undef, I,J,N)
    Ic00 = BitArray{3}(undef, I,J,N)

    K_supply  = 0.
    L_supply  = 0.
    KL_supply = 0.
    Vamin     = 0.
    Vbmin     = 1e-8

    #----------------------------------------------------------------
    # Iterate on KL to find steady state
    #----------------------------------------------------------------
    for ii = 1 : maxit_KL
	    # Derive aggregates, given KL
	    w   = (1 - aalpha) * (KL ^ aalpha)
 	    r_a = aalpha * (KL ^ (aalpha - 1)) - ddelta

	    # Store current value function
	    Vn = V_0

	    #----------------------------------------------------------------
	    # Solve HJB
	    #----------------------------------------------------------------
        @time for nn = 1 : maxit_HJB
            perm_c =  ddeath * pam # NOTE: perm?
            c0_c   = ((1-xxi) - tau_I) * w

            for i=1:I, j=1:J, n=1:N
                c0 = c0_c * real(y[n]) + b[i] * (r_b_vec[i] + perm_c) + trans

                # ---- Liquid Assets, Forward + Backward Difference ----
                VaF = (j==J) ? 0.0 : max((Vn[i,j+1,n] - Vn[i,j,n]) / (a[j+1]-a[j]), Vamin)
                VaB = (j==1) ? 0.0 : max((Vn[i,j,n] - Vn[i,j-1,n]) / (a[j]-a[j-1]), Vamin)

                # ---- Illiquid Assets, Forward + Backward Difference ----
                VbF = (i==I) ? 0.0 : max((Vn[i+1,j,n] - Vn[i,j,n]) / (b[i+1]-b[i]), Vbmin)
                VbB = (i==1) ? 0.0 : max((Vn[i,j,n] - Vn[i-1,j,n]) / (b[i]-b[i-1]), Vbmin)

                # Decisions conditional on a particular direction of derivative
                cF = (i==I) ? 0.0 : VbF ^ (-1 / ggamma)
                cB = (i==1) ? c0  : VbB ^ (-1 / ggamma)

                sF = (i==I) ? 0.0 : c0 - cF
                sB = (i==1) ? 0.0 : c0 - cB

                if (i==1 && j > 1) VbB = util(cB) end

                #----------------------------------------------------------------
                # Consumption  & Savings Decision
                #----------------------------------------------------------------
                # Which direction to use
                Hc0 = util(c0)
                HcF = (i==I) ? -1e12 : util(cF) + VbF * sF
                HcB = util(cB) + VbB * sB

                validF = (sF > 0)
                validB = (sB < 0)

                IcF[i,j,n] = validF * max(!validB, (HcF >= HcB)) * (HcF >= Hc0)
                IcB[i,j,n] = validB * max(!validF, (HcB >= HcF)) * (HcB >= Hc0)
                Ic0[i,j,n] = 1 - IcF[i,j,n] - IcB[i,j,n]

                c[i,j,n] = IcF[i,j,n] * cF + IcB[i,j,n] * cB + Ic0[i,j,n] * c0
                s[i,j,n] = IcF[i,j,n] * sF + IcB[i,j,n] * sB

                #----------------------------------------------------------------
                # Deposit Decision
                #----------------------------------------------------------------
                dFB = (i == 1 || j == J) ? 0.0 : deposit(VaF, VbB, max(a[j], a_lb))
                dBF = (i == I || j == 1) ? 0.0 : deposit(VaB, VbF, max(a[j], a_lb))
                dBB = (j == 1)           ? 0.0 : deposit(VaB, VbB, max(a[j], a_lb))

                HdFB = (i == 1 || j == J) ? -1e12 : VaF * dFB - VbB * (dFB + cost(dFB, a[j]))
                validFB = (dFB > 0) * (HdFB > 0)

                HdBB = (j == 1) ? -1.0e12 : VaB * dBB - VbB * (dBB + cost(dBB, a[j]))
                validBB = (dBB > -cost(dBB, a[j])) .* (dBB <= 0) * (HdBB > 0)

                HdBF = (i == I || j == 1) ? -1e12 : VaB * dBF - VbF * (dBF + cost(dBF, a[j]))
                validBF = (dBF <= -cost(dBF, a[j])) * (HdBF > 0)

                IcFB[i,j,n] = validFB * max(!validBF, HdFB >= HdBF) * max(!validBB, HdFB >= HdBB)
                IcBF[i,j,n] = max(!validFB, HdBF >= HdFB) * validBF * max(!validBB, HdBF >= HdBB)
                IcBB[i,j,n] = max(!validFB, HdBB >= HdFB) * max(!validBF, HdBB >= HdBF) * validBB
                Ic00[i,j,n] = !validFB * !validBF * !validBB

                d[i,j,n] = IcFB[i,j,n] * dFB + IcBF[i,j,n] * dBF + dBB * IcBB[i,j,n]
            end
            u = util.(c)

            # Interpolate
            d_g = reshape(interp_decision * vec(d), I_g, J_g, N)
            s_g = reshape(interp_decision * vec(s), I_g, J_g, N)
            c_g = reshape(interp_decision * vec(c), I_g, J_g, N)

            aa, bb, aau, bbu = transition(ddeath, pam, xxi, w, chi0, chi1, chi2, a_lb, r_a,
                                          vec(y), d, d_g, s, s_g, a, a_g, b, b_g)
            A = aa + bb

            #------------------------------------------------------------
            # Update value function
            #------------------------------------------------------------
            Vn1 = Array{Float64}(undef, I,J,N)
            for kk = 1:N
                Ak          = A[1+(kk-1)*(I*J):kk*(I*J), 1+(kk-1)*(I*J):kk*(I*J)]
                Bk          = (1 + Delta*(rrho + ddeath) -
                               Delta*lambda[kk,kk]) * my_speye(I*J) - Delta*Ak
                uk_stacked  = reshape(u[:,:,kk], I * J, 1)
                Vk_stacked  = reshape(Vn[:,:,kk], I * J, 1)
                indx_k      = 1:N .!= kk
                Vkp_stacked = sum(broadcast(*,  lambda[kk, indx_k]',
                                          reshape(Vn[:,:,indx_k],I*J,N-1)), dims=2)
                qk          = Delta * uk_stacked + Vk_stacked + Delta * Vkp_stacked
                Vn1k_stacked = Bk \ qk
                Vn1[:,:,kk]  = reshape(Vn1k_stacked, I, J, 1)
            end

            # Howard improvement step
            if nn >= start_HIS
                for jj = 1:maxit_HIS
                    Vn2 = Array{Float64}(undef, I, J, N)
                    for kk = 1:N
                        uk_stacked = reshape(u[:,:,kk],I*J,1)
                        Vk_stacked = reshape(Vn1[:,:,kk],I*J,1)
                        Ak = A[1+(kk-1)*(I*J):kk*(I*J), 1+(kk-1)*(I*J):kk*(I*J)]
                        Bk = (1 + Delta*(rrho + ddeath) -
                              Delta*lambda[kk,kk]) * my_speye(I*J) - Delta*Ak
                        indx_k         = 1:N .!= kk
                        Vkp_stacked = sum(broadcast(*,  lambda[kk, indx_k]',
                                          reshape(Vn1[:,:,indx_k],I*J,N-1)), dims=2)
                        qk             = Delta*uk_stacked + Vk_stacked + Delta*Vkp_stacked
                        Vn2k_stacked = Bk \ qk
                        Vn2[:,:,kk] = reshape(Vn2k_stacked,I,J,1)
                    end
                    VHIS_delta = Vn2 - Vn1
                    Vn1  = Vn2
                    dist = maximum(abs.(VHIS_delta))

                    if dist < crit_HIS
                        break
                    end
                end
            end

            # Check for convergence
            V_Delta = Vn1 - Vn
            Vn      = Vn1

            dist = maximum(abs.(V_Delta))
            @show dist

            if dist < crit_HJB
                 println("Value Function Converged, Iteration = ", nn)
                break
            end
        end

        # Store value function
        V_0 = Vn

        #----------------------------------------------------------------
        # Compute stationary distribution and update aggregate KL
        #----------------------------------------------------------------

        # Find new stationary distribution associated with decision rules
        A   = aau + bbu
        λ0  = lambda - diagm(0 => diag(lambda))  # transition matrix with diagonal killed
        λ0p = λ0'

        gg_tilde = dab_g_tilde_mat * gg
        gg1      = Array{Float64,2}(undef, I_g * J_g, N)
        g        = Array{Float64,3}(undef, I_g, J_g, N)

        K_supply  = 0.
        L_supply  = 0.
        KL_supply = 0.

        # Iterate
        for nn = 1:maxit_KFE

            gg_tilde = dab_g_tilde_mat * gg
            gg1      = Array{Float64,2}(undef, I_g * J_g, N)

            for kk = 1:N
                Ak = A[1+(kk-1)*(I_g*J_g):kk*(I_g*J_g),
                       1+(kk-1)*(I_g*J_g):kk*(I_g*J_g)]

                death_inflow = zeros(I_g, J_g)
                death_inflow[b .== 0, 1] .= sum(gg_tilde[1+(kk-1)*(I_g*J_g):kk*(I_g*J_g)])

                death_inflow = reshape(death_inflow, I_g*J_g, 1)

                gk_sum = sum(repeat(λ0p[kk,:]', I_g*J_g,1) .* reshape(gg_tilde,I_g*J_g,N), dims=2)
                gg1[:,kk] = (my_speye(I_g*J_g) - Delta_KFE * Ak' - Delta_KFE *
                             (lambda[kk,kk] - ddeath) *
                             my_speye(I_g*J_g))\(gg_tilde[1+(kk-1)*(I_g*J_g):kk*(I_g*J_g)] +
                                                 Delta_KFE*gk_sum + Delta_KFE*ddeath*death_inflow)
            end

            gg1     = reshape(gg1, I_g*J_g*N,1)
            gg1     = gg1 ./ sum(gg1)
            gg1     = dab_g_tilde_mat \ gg1

            dist = maximum(abs.(gg1-gg))

            if mod(nn,100) == 0
                println("Iteration ", nn, ", distance is ", dist)
            end
            if dist < crit_KFE
                gg = gg1
                g  = reshape(gg, I_g, J_g, N)

                println("Distribution Converged, Iteration = ", nn)
                break
            end
            gg = gg1
            g  = reshape(gg, I_g, J_g, N)
        end

        #------------------------------------------------------------
        # Update guess of KL
        #------------------------------------------------------------
        # Capital supply
        if K_liquid == 1
            K_supply = sum(g .* (a_g_grid + b_g_grid) .* dab_tilde_grid)
        else
            K_supply = sum(g .* a_g_grid .* dab_g_tilde_grid)
        end

        # Labor supply
        L_supply = sum(y_g_grid .* g .* dab_g_tilde_grid)

        # Capital-Labor ratio
        KL_supply = K_supply / L_supply

        # Check for convergence and update
        gap = KL - KL_supply
        println("The current gap is ", gap)

        if abs(gap) > crit_KL
            KL = relax_KL * KL_supply + (1 - relax_KL) * KL
        else
            println("I have found the steady state, Iteration = ", ii)
            break
        end

    end

    #compute_savings()
    ccu                  = kron(lambda, my_speye(I_g * J_g))
    A                    = aau + bbu + ccu
    dab_small            = reshape(dab_g_tilde_grid[:,:,1], I_g*J_g, 1)
    loc                  = findall(!iszero, b .== 0)
    dab_small            = dab_small ./ dab_small[loc] * ddeath
    dab_small[loc]      .= 0.0
    death_process        = -ddeath * my_speye(I_g*J_g)
    death_process[loc,:] = vec(dab_small)
    death_process        = kron(my_speye(N), death_process)

    g_new = A' * gg_tilde
    g_new = dab_g_tilde_mat \ g_new
    g_new = g_new + death_process * gg
    g_new = reshape(g_new ,I_g, J_g, N)
    g_new_a = sum(sum(g_new, dims=1), dims=3)
    save_a = dot(g_new_a, a_g)
    #g_new_b = sum(sum(g_new, dims=2), dims=3)
    #save_b = dot(vec(g_new_b), b_g)
    #compute_savings()

    # Rename variables in steady state
    m[:KL_SS]           = real(KL)
    m[:r_a_SS]          = aalpha * (m[:KL_SS].value ^ (aalpha - 1)) - ddelta
    m[:w_SS]            = (1 - aalpha) * (m[:KL_SS].value ^ aalpha)
    m[:K_SS]            = real(m[:KL_SS].value * L_supply)
    m[:u_SS]            = u
    m[:c_SS]            = c
    m[:d_SS]            = d
    m[:V_SS]            = Vn
    m[:g_SS]            = g
    m[:dab]             = dab_tilde_grid
    m[:dab_g]           = dab_g_tilde_grid
    m[:C_SS]            = sum(g .* c_g .* dab_g_tilde_grid)
    m[:C_Var_SS]        = sum(g .* log.(c_g).^2 .* dab_g_tilde_grid) -
                                 sum(g .* log.(c_g) .* dab_g_tilde_grid)^2
    m[:I_SS]            = save_a
    m[:B_SS]            = sum(g .* b_g_grid .* dab_g_tilde_grid)
    m[:Y_SS]            = real((K_supply ^ aalpha) * (L_supply ^ (1 - aalpha)))
    m[:n_SS]            = real(L_supply)
    m[:earn_SS]         = log.((1-tau_I) * w * y_g_grid + b_g_grid .*
                               (r_b_g_grid .+ ddeath*pam) .+ trans .+ a_g_grid .*
                               (r_a .+ ddeath*pam))
    m[:earn_Var_SS]     = sum(g .* m[:earn_SS].value .^ 2 .* dab_g_tilde_grid) -
                                sum(g .* m[:earn_SS].value .* dab_g_tilde_grid)^2
    m <= Setting(:IcF_SS,  IcF)
    m <= Setting(:IcB_SS,  IcB)
    m <= Setting(:Ic0_SS,  Ic0)
    m <= Setting(:IcFB_SS, IcFB)
    m <= Setting(:IcBF_SS, IcBF)
    m <= Setting(:IcBB_SS, IcBB)
    m <= Setting(:Ic00_SS, Ic00)

    a_g_0pos = get_setting(m, :a_g_0pos)
    b_g_0pos = get_setting(m, :b_g_0pos)

    ###
    # Consumption by hand-to-mouth status
    ###
    WHTM_indicator      = zeros(I_g, J_g, N)
    WHTM_indicator[b_g_0pos:b_g_0pos+1, a_g_0pos+2:end,:] .= 1.0

    m[:WHTM_indicator] = WHTM_indicator
    m[:WHTM_SS]        = sum(g[:] .* WHTM_indicator[:] .* dab_g_tilde_grid[:])
    m[:C_WHTM_SS]      = sum(WHTM_indicator[:] .* c_g[:] .* g[:] .* dab_g_tilde_grid[:])

    PHTM_indicator     = zeros(I_g, J_g, N)
    PHTM_indicator[b_g_0pos:b_g_0pos+1,a_g_0pos:a_g_0pos+2:end,:] .= 1.0
    m[:PHTM_indicator] = PHTM_indicator
    m[:PHTM_SS]        = sum(g[:] .* PHTM_indicator[:] .* dab_g_tilde_grid[:])
    m[:C_PHTM_SS]      = sum(c_g[:] .* g[:] .* PHTM_indicator[:] .* dab_g_tilde_grid[:])

    NHTM_indicator      = zeros(I_g,J_g,N)
    NHTM_indicator[(b_g_0pos+3):end,2:end,:] .= 1.0

    m[:NHTM_indicator] = NHTM_indicator
    m[:NHTM_SS]        = sum(g[:] .* NHTM_indicator[:] .* dab_g_tilde_grid[:])
    m[:C_NHTM_SS]      = sum(c_g[:] .* g[:] .* NHTM_indicator[:] .* dab_g_tilde_grid[:]) /
        m[:NHTM_SS].value
    m[:r_a_grid]       = repeat([r_a], I_g, J_g, N)

    n_v = get_setting(m, :n_v)
    n_g = get_setting(m, :n_g)
    n_p = get_setting(m, :n_p)
    n_Z = get_setting(m, :n_Z)

    # Collect variables into vector
    vars_SS                  = zeros(n_v + n_g + n_p + n_Z, 1)
    vars_SS[1:n_v]           = reshape(m[:V_SS].value, I*J*N, 1)

    gg_SS                    = reshape(m[:g_SS].value, I_g*J_g*N, 1)
    m[:gg_SS]                = gg_SS

    vars_SS[n_v+1:n_v+n_g,1] = gg_SS[1:I_g*J_g*N-1]
    vars_SS[n_v+n_g+1,1]     = m[:K_SS].value
    vars_SS[n_v+n_g+2,1]     = m[:r_b_SS].value

    if aggregate_variables == 1
        vars_SS[n_v+n_g+3,1]         = m[:Y_SS].value
        vars_SS[n_v+n_g+4,1]         = m[:C_SS].value
    elseif distributional_variables == 1
        vars_SS[n_v+n_g+3,1]        = m[:C_Var_SS].value
        vars_SS[n_v+n_g+4,1]        = m[:earn_Var_SS].value
    elseif distributional_variables_1 == 1
        vars_SS[n_v+n_g+3,1]        = m[:C_WHTM_SS].value
        vars_SS[n_v+n_g+4,1]        = m[:C_PHTM_SS]
        #vars_SS[n_v+n_g+3,1]        = C_NHTM_SS
    end
    vars_SS[n_v+n_g+n_p+1,1] = 0.0

    # Compute illiquid wedge
    m[:illiquid_wedge] = m[:r_a_SS].value - m[:r_b_SS].value

    # Save SS variables
    m[:vars_SS] = vars_SS

    return m
end

"""
```
function adjust_parameters!(m::TwoAssetHANK)
```
Function adjusts model parameters depending on what settings are provided.
"""
function adjust_parameters!(m::TwoAssetHANK)
    if get_setting(m, :para_new) == 0
        m[:rrho]   = 0.014526858925819 # discount factor (quarterly)
        m[:ddeath] = 1. / (4. * 45)    # death rate (quarterly)
        # Technology
        m[:aalpha] = 0.4     # capital share
        m[:ddelta] = 0.02    # depreciation rate (quarterly)

        # Portfolio adjustment cost function
        m[:chi0] = 0.070046
        m[:chi1] = 0.535791
        m[:chi2] = 0.971050
        m[:a_lb] = 0.2407

        # Tax code
        m[:trans] = 0.430147477922432 # transfer to households
        m[:tau_I] = 0.25              # tax rate on wage income

        # Liquid return
        m[:borrwedge_SS] = 0.019875
        m[:r_b_borr_SS]  = m[:r_b_SS] + m[:borrwedge_SS]

        # Stochastic Properties of aggregate TFP shock
        if !get_setting(m, :permanent)
            m[:nnu_aggZ]   = 0.05  # persistence parameter in Ornstein-Uhlenbeck process
            m[:sigma_aggZ] = 0.007 # volatility parameter in Ornstein-Uhlenbeck process
        else
            m[:nnu_aggZ]   = 2.15 # persistence parameter in Ornstein-Uhlenbeck process
            m[:sigma_aggZ] = 0.04 # volatility parameter in Ornstein-Uhlenbeck process
        end
    else

        m[:rrho]   = 0.014495841360662
        m[:ddelta] = 0.075 / 4
        m[:a_lb]   = 0.02
        m[:trans] = 0.1
        m[:tau_I] = 0.3
        m[:xxi]   = 0.0

        # Stochastic properties of aggregate TFP shock
        if !get_setting(m, :permanent)
            m[:nnu_aggZ]    = 0.25
            m[:ssigma_aggZ] = 0.007 / (1.0 - m[:aalpha])
        else
            m[:nnu_aggZ]    = 1.265
            m[:ssigma_aggZ] = 0.0282
        end
    end

    # Household effects of TFP shock
    if get_setting(m, :asymm)
        if get_setting(m, :y_size) == 2
            m[:kappa] = 10.0
        elseif get_setting(m, :y_size) == 33
            m[:kappa] = 2.0
        end
    else
        m[:kappa] = 0.0
    end

    # Necessary settings for permanent shock
    if get_setting(m, :permanent)
        m <= Setting(:reduceDist_hor, 400)
        m <= Setting(:r_b_fix, true)
    end
end

function model_settings!(m::TwoAssetHANK)
    default_settings!(m)

    # Productivity shock properties
    m <= Setting(:permanent, false,
                 "Permanent or transitory shock. For permanent, use reduceDist_hor=400, r_b_fix=1.")
    m <= Setting(:asymm,     false,
                 "Does effect of shock vary according to position in the wealth distribution?")

    # Solution approach
    m <= Setting(:reduceDistribution, true, "Reduction of g using observability matrix")
    m <= Setting(:reduceV,            true, "Reduction of v using splines")
    m <= Setting(:reduceDist_hor,     100,
                 "Maximal horizon for observability matrix")
    m <= Setting(:krylov_dim, 100, "Maximal horizon for observability matrix") # same as above


    # Liquid asset suply
    m <= Setting(:r_b_fix,  true, "Adjustment procedure 1: r_b is fixed")
    m <= Setting(:r_b_phi,  false,
                 "Adjustment procedure 2: r_b is supplied elastically, with elasticity phi")
    m <= Setting(:B_fix,    false, "The total amount of the liquid asset is fixed at B")
    m <= Setting(:K_liquid, false, "K_liquid")

    # Model size
    m <= Setting(:y_size, 30, "Model size: number of income states. Supported sizes: 2, 30, 33")

    # Distributional variables
    m <= Setting(:aggregate_variables,        true, "Aggregate variables")
    m <= Setting(:distributional_variables,   false, "Distributional_variables")
    m <= Setting(:distributional_variables_1, false, "Distributional_variables_1")

    # Parameterization
    m <= Setting(:para_new, 1, "Parameterization")

    # Grids
    m <= Setting(:agrid_new, 0, "agrid_new = 0")
    m <= Setting(:bgrid_new, 0, "bgrid_new = 0")
    m <= Setting(:ygrid_new, 0, "ygrid_new = 0")

    # Create income grid
    lambda, y, y_mean, y_dist = create_y_grid(get_setting(m, :y_size),
                                              get_setting(m, :ygrid_new))
    m <= Setting(:lambda, lambda, "Transition probabilities")
    m <= Setting(:y, y, "Grid of income shocks")
    m <= Setting(:y_mean, y_mean, "Mean grid of income shocks")
    m <= Setting(:y_dist, y_dist, "Stationary distribution of income shocks")
    m <= Setting(:N, length(y), "N")

    # Create illiquid asset grid
    m <= Setting(:J,    40, "Number of illiquid grid points (a)")
    m <= Setting(:J_g,  get_setting(m, :agrid_new) == 0 ? 45  : 40, "")
    m <= Setting(:amax, get_setting(m, :agrid_new) == 0 ? 200 : 2000,
                 "Maximum illiquid asset grid value")
    m <= Setting(:amin, 0, "Minimum illiquid asset grid value")
    m <= Setting(:a_,   0, "Lower bound on illiquid assets")

    a, a_g, a_g_0pos = create_a_grid(get_setting(m, :agrid_new), get_setting(m, :J),
                                     get_setting(m, :J_g), get_setting(m, :amin),
                                     get_setting(m, :amax))
    m <= Setting(:a,        a, "a grid")
    m <= Setting(:a_g,      a_g, "a_g grid")
    m <= Setting(:a_g_0pos, a_g_0pos[1], "Points in illiquid asset grid equal to 0")

    # Create liquid asset grid
    m <= Setting(:I,   50, "Number of liquid grid points")
    m <= Setting(:I_g, 50, "I_g")

    I_neg, bmin, bmax, b, b_g, b_g_0pos = create_b_grid(get_setting(m, :agrid_new),
                                                        get_setting(m, :I), get_setting(m, :I_g))
    m <= Setting(:I_neg,    I_neg, "Number of liquid grid points on negative part of real line")
    m <= Setting(:I_pos,    get_setting(m, :I) - I_neg,
                 "Number of liquid grid points on positive part of real line")
    m <= Setting(:bmin,     bmin, "Lower bound on liquid assets")
    m <= Setting(:bmax,     bmax, "Upper bound on liquid assets")
    m <= Setting(:b,        b, "b")
    m <= Setting(:b_g,      b_g, "b_g")
    m <= Setting(:b_g_0pos, b_g_0pos[1], "Points in liquid asset grid equal to 0")

    m <= Setting(:interp_decision, kron(my_speye(get_setting(m, :N)), interp(b_g, a_g, b, a)),
                                        "interp_decision")

    m <= Setting(:n_v,      get_setting(m, :I) * get_setting(m, :J) * get_setting(m, :N),
                 "Number of jump variables (value function entries)")
    m <= Setting(:n_g,      get_setting(m, :I_g) * get_setting(m, :J_g) * get_setting(m, :N) - 1,
                 "Number of endogenous state variables (distribution)")
    m <= Setting(:n_p,     (get_setting(m, :aggregate_variables) == 1) ||
                           (get_setting(m, :distributional_variables) == 1) ? 4 : 2, "n_p")
    m <= Setting(:n_Z,      1, "Number of aggregate shocks")
    m <= Setting(:nVars,    get_setting(m, :n_v) + get_setting(m, :n_g) + get_setting(m, :n_p) +
                            get_setting(m, :n_Z), "nVars")
    m <= Setting(:nEErrors, get_setting(m, :n_v), "nEErrors")
    m <= Setting(:n_v_full, get_setting(m, :n_v), "n_v_Full")
    m <= Setting(:n_g_full, get_setting(m, :n_g), "n_g_Full")

    # Our syntax
    m <= Setting(:n_jump_vars,  get_setting(m, :I) * get_setting(m, :J) * get_setting(m, :N),
                 "Number of jump variables (value function entries)")
    m <= Setting(:n_state_vars, get_setting(m, :I_g)*get_setting(m, :J_g)*get_setting(m, :N)-1,
                 "Number of endogenous state variables (distribution)")
    m <= Setting(:n_state_vars_unreduce, 1, "Number of aggregate shocks")

    # Steady state approximation
    m <= Setting(:maxit_HJB, 100,  "Max number of iterations for HJB")
    m <= Setting(:crit_HJB,  1e-6, "Tolerance for HJB error")
    m <= Setting(:Delta,     1e6,  "Multiplier for implicit upwind scheme of HJB")

    m <= Setting(:maxit_HIS, 10,   "Max allowable number of Howard improvement steps")
    m <= Setting(:start_HIS, 2,    "When in HJB loop should Howard improvement steps begin?")
    m <= Setting(:crit_HIS,  1e-5, "Critical value HIS")

    m <= Setting(:maxit_KFE, 1000, "Max allowable number of KFE iterations")
    m <= Setting(:crit_KFE,  1e-7, "Tolerance for KFE error")
    m <= Setting(:Delta_KFE, 1000, "Step size KFE")

    m <= Setting(:maxit_KL,  1000, "Max number of iterations over capital-labor ratio")
    m <= Setting(:crit_KL,   1e-5, "Tolerance for KL error")
    m <= Setting(:relax_KL,  0.1,  "Updating speed (must lie within [0,1])")

    # Approximation Parameters
    m <= Setting(:KL_0, 43.8800, "Initial guess of capital to labor ratio")

    # TODO: Inserted from old code; not certain necessary
    #=

    m <= Setting(:n_static_relations, 5,
                 "Number of static relations: bond-market clearing, labor market clearing, consumption, output, total assets")
    m <= Setting(:n_vars, get_setting(m, :n_jump_vars) + get_setting(m, :n_state_vars)
                 + get_setting(m, :n_static_relations),
                 "Number of variables, total")

    # Reduction
    m <= Setting(:n_knots, 12, "Number of knot points")
    m <= Setting(:c_power, 1, "Amount of bending of knot point locations to make them nonuniform")
    m <= Setting(:n_post, length(get_setting(m, :zz)[1,:]),
                 "Number of dimensions that need to be approximated by spline basis")
    m <= Setting(:n_prior, 1,
                 "Number of dimensions approximated by spline basis that
                  were not used to compute the basis matrix")

    knots = collect(range(get_setting(m, :amin), stop=get_setting(m, :amax),
                          length=get_setting(m, :n_knots)-1))
    knots = (get_setting(m, :amax) - get_setting(m, :amin)) /
        (2^get_setting(m, :c_power)-1) * ((knots .- get_setting(m, :amin))
             / (get_setting(m, :amax) - get_setting(m, :amin)) .+ 1) .^
             get_setting(m, :c_power) .+ get_setting(m, :amin) .- (get_setting(m, :amax) -
             get_setting(m, :amin)) / (2 ^ get_setting(m, :c_power) - 1)
    m <= Setting(:knots_dict, Dict(1 => knots),
                 "Location of knot points for each dimension for value function reduction")
    m <= Setting(:krylov_dim, 20, "Krylov reduction dimension")
    m <= Setting(:reduce_state_vars, true, "Reduce state variables or not") #turned off reduction
    m <= Setting(:reduce_v, true, "Reduce value function or not") #turned off reduction
    m <= Setting(:spline_grid, get_setting(m, :a), "Grid of knot points for spline basis")
=#
    # Sampling method
    m <= Setting(:sampling_method, :SMC, "Set sampling method to SMC")

    # SMC Settings
    m <= Setting(:n_particles, 3000)
    m <= Setting(:n_Φ, 200)
    m <= Setting(:λ, 2.0)
    m <= Setting(:n_smc_blocks, 1)
    m <= Setting(:use_parallel_workers, true)
    m <= Setting(:step_size_smc, 0.5)
    m <= Setting(:n_mh_steps_smc, 1)
    m <= Setting(:resampler_smc, :polyalgo)
    m <= Setting(:target_accept, 0.25)
    m <= Setting(:mixture_proportion, .9)
    m <= Setting(:tempering_target, 0.95)
    m <= Setting(:resampling_threshold, .5)
    m <= Setting(:use_fixed_schedule, true)
    m <= Setting(:smc_iteration, 0)

    # Track lag
    m <= Setting(:track_lag, false, "Add first lag when constructing measurement equation")

    # Forecast
    m <= Setting(:use_population_forecast, true,
                 "Whether to use population forecasts as data")
    m <= Setting(:forecast_zlb_value, 0.13,
        "Value of the zero lower bound in forecast periods, if we choose to enforce it")

    # Simulating states
    m <= Setting(:state_simulation_freq, 2,
                 "How many states you want to simulate between states + 1")
end
