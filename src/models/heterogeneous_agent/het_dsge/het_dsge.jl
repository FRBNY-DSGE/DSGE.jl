import DataStructures: OrderedDict

"""
```
HetDSGE{T} <: AbstractModel{T}
```

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

* `observables::OrderedDict{Symbol,Int}`: Maps each observable to a row in the
  model's measurement equation matrices.

#### Model Specifications and Settings

* `spec::String`: The model specification identifier, \"an_schorfheide\", cached
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

"""
mutable struct HetDSGE{T} <: AbstractModel{T}
    parameters::ParameterVector{T}                         # vector of all time-invariant model parameters
    steady_state::ParameterVector{T}                       # model steady-state values

    # Temporary to get it to work. Need to
    # figure out a more flexible way to define
    # "grids" that are not necessarily quadrature
    # grids within the model
    grids::OrderedDict{Symbol,Union{Grid, Array}}

    keys::OrderedDict{Symbol,Int}                          # human-readable names for all the model
                                                           # parameters and steady-states

    state_variables::Vector{Symbol}                        # Vector of symbols of the state variables
    jump_variables::Vector{Symbol}                         # Vector of symbols of the jump variables
    normalized_model_states::Vector{Symbol}                # All of the distributional model
                                                           # state variables that need to be normalized
    endogenous_states_unnormalized::OrderedDict{Symbol,UnitRange} # Vector of unnormalized
                                                           # ranges of indices
    endogenous_states::OrderedDict{Symbol,UnitRange}       # Vector of ranges corresponding
                                                           # to normalized (post Klein solution) indices
    exogenous_shocks::OrderedDict{Symbol,Int}              #
    expected_shocks::OrderedDict{Symbol,Int}               #
    equilibrium_conditions::OrderedDict{Symbol,UnitRange}  #
    observables::OrderedDict{Symbol,Int}                   #

    spec::String                                           # Model specification number (eg "m990")
    subspec::String                                        # Model subspecification (eg "ss0")
    settings::Dict{Symbol,Setting}                         # Settings/flags for computation
    test_settings::Dict{Symbol,Setting}                    # Settings/flags for testing mode
    rng::MersenneTwister                                   # Random number generator
    testing::Bool                                          # Whether we are in testing mode or not

    observable_mappings::OrderedDict{Symbol, Observable}
end

description(m::HetDSGE) = "HetDSGE, $(m.subspec)"

"""
`init_model_indices!(m::HetDSGE)`

Arguments:
`m:: HetDSGE`: a model object

Description:
Initializes indices for all of `m`'s states, shocks, and equilibrium conditions.
"""
function init_model_indices!(m::HetDSGE)
    # Endogenous states
    endogenous_states = collect([:l′_t1, :μ′_t1, :k′_t, :R′_t1, :i′_t1,
                                 :π′_t1,:π′_t2,:π′_t3,:y′_t1,:y′_t2,:y′_t3, :y′_t4,
                                 :z′_t1,:z′_t2,:z′_t3,:w′_t1,:I′_t1,:BP,:GP,:ZP,:MUP,:LAMWP,
                                 :LAMFP,:MONP,:l′_t,:μ′_t,:R′_t,:i′_t,:t′_t,:w′_t,:L′_t,:π′_t,
                                 :wageinflation′_t,:mu′_t,:y′_t,:I′_t,:mc′_t,:Q′_t,:capreturn′_t,
                                 :l_t1L,:μ_t1,:k_t,:R_t1,:i_t1,:π_t1,:π_t2,:π_t3,:y_t1,:y_t2,
                                 :y_t3,:y_t4,:z_t1,:z_t2,:z_t3,:w_t1,:I_t1,:B,:G,:Z,:MU,:LAMW,
                                 :LAMF,:MON,:ELL,:μ_t,:R_t,:i_t,
                                 :t_t,:w_t,:L_t,:π_t,:wageinflation_t,:mu_t,:y_t,:I_t,:mc_t,:Q_t,:capreturn_t])

    # Exogenous shocks
    exogenous_shocks = collect([])

    # Equilibrium conditions
    equilibrium_conditions = collect([:eq_euler,:eq_kolmogorov_fwd,:eq_lag_ell,:eq_lag_wealth,
                                      :eq_market_clearing,:eq_lambda, :eq_transfers,
                                      :eq_investment,:eq_tobin_q,:eq_capital_accumulation,
                                      :eq_wage_phillips,:eq_price_phillips,:eq_marginal_cost,
                                      :eq_gdp,:eq_optimal_kl,:eq_taylor, :eq_fisher,
                                      :eq_nominal_wage_inflation,:LR,:LI,:LPI,:L2PI,:L3PI,:LY,:L2Y,
                                      :L3Y,:L4Y,:LZ,:L2Z,:L3Z,:LW,:LX,:F33,:F34,:F35,:F36,:F37,
                                      :F38, :F39])

    # Observables
    observables = keys(m.observable_mappings)

    ########################################################################################
    # Setting indices of endogenous_states and equilibrium conditions manually for now

    setup_indices!(m)
    endo = m.endogenous_states_unnormalized
    eqcond = equilibrium_conditions

    ########################################################################################

    m.normalized_model_states = [:μ′_t]
    m.state_variables = m.endogenous_states.keys[get_setting(m, :state_indices)]
    m.jump_variables = m.endogenous_states.keys[get_setting(m, :jump_indices)]

    for (i,k) in enumerate(exogenous_shocks);            m.exogenous_shocks[k]            = i end
    for (i,k) in enumerate(observables);                 m.observables[k]                 = i end
end

function HetDSGE(subspec::String="ss0";
                   custom_settings::Dict{Symbol, Setting} = Dict{Symbol, Setting}(),
                   testing = false)

    # Model-specific specifications
    spec               = "het_dsge"
    subspec            = subspec
    settings           = Dict{Symbol,Setting}()
    test_settings      = Dict{Symbol,Setting}()
    rng                = MersenneTwister(0)

    # initialize empty model
    m = HetDSGE{Float64}(
            # model parameters and steady state values
            Vector{AbstractParameter{Float64}}(), Vector{Float64}(),
            # grids and keys
            OrderedDict{Symbol,Union{Grid, Array}}(), OrderedDict{Symbol,Int}(),

            # normalized_model_states, state_inds, jump_inds
            Vector{Symbol}(), Vector{Symbol}(), Vector{Symbol}(),

            # model indices
            # endogenous states unnormalized, endogenous states normalized
            OrderedDict{Symbol,UnitRange}(), OrderedDict{Symbol,UnitRange}(),
            OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(),
            OrderedDict{Symbol,UnitRange}(), OrderedDict{Symbol,Int}(),

            spec,
            subspec,
            settings,
            test_settings,
            rng,
            testing,
            OrderedDict{Symbol,Observable}())

    # Set settings
    model_settings!(m)
    # default_test_settings!(m)
    for custom_setting in values(custom_settings)
        m <= custom_setting
    end

    # # Set observable transformations
    # init_observable_mappings!(m)

    # Initialize parameters
    init_parameters!(m)

    # Initialize aggregate steady state parameters (necessary for grid construction)
    aggregate_steadystate!(m)

    # Initialize grids
    init_grids!(m)

    # Initialize model indices
    init_model_indices!(m)

    # # Solve for the steady state
    # steadystate!(m)

    # # So that the indices of m.endogenous_states reflect the normalization
    # normalize_model_state_indices!(m)

    return m
end

"""
```
init_parameters!(m::HetDSGE)
```

Initializes the model's parameters, as well as empty values for the steady-state
parameters (in preparation for `steadystate!(m)` being called to initialize
those).
"""
function init_parameters!(m::HetDSGE)

    # Initialize parameters
    m <= parameter(:r, 0.01, fixed = true,
                   description = "r: Steady-state real interest rate.", tex_label = "r")
    m <= parameter(:α, 0.3, fixed = true,
                   description = "α: Capital share", tex_label = "\\alpha")
    m <= parameter(:H, 1.0, fixed = true,
                   description = "H: Aggregate hours worked", tex_label = "H")
    m <= parameter(:δ, 0.03, fixed = true,
                   description = "δ: Depreciation rate", tex_label = "\\delta")
    m <= parameter(:μ_sp, 0.0, fixed = true, description = "μ_sp: The trend in the skill process",
                   tex_label = "\\mu_{sp}")
    m <= parameter(:ρ_sp, 0.7, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), BetaAlt(0.5, 0.2), fixed=false,
                   description="ρ_sp: AR(1) coefficient in the skill process.",
                   tex_label="\\rho_{sp}")
    m <= parameter(:σ_sp, 0.01, (1e-8, 5.), (1e-8, 5.), Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_sp: The standard deviation of the skill process.",
                   tex_label="\\sigma_{sp}")
    m <= parameter(:γ, 0.0, description = "γ: TFP growth")
    #GoverY = 0.01
    m <= parameter(:g, 1/(1-0.01), description = "g: Steady-state government spending/gdp")
    m <= parameter(:η, 0.1, description = "η: Borrowing constraint (normalized by TFP)")

    m <= parameter(:ρB, 0.5, description = "# persistence of discount factor shock")
    m <= parameter(:ρG, 0.5, description = "# persistence of govt spending shock")
    m <= parameter(:ρZ, 0.5, description = "# persistence of tfp growth shock")
    m <= parameter(:ρμ, 0.5, description = " # persistence of investment shock")
    m <= parameter(:ρlamw, 0.5, description = "# persistence of wage mkup shock")
    m <= parameter(:ρlamf, 0.5, description = " # persistence of price mkup shock")
    m <= parameter(:ρmon, 0.5, description = " # persistence of mon policy shock")

 m <= parameter(:spp, 4., description = "# second derivative of investment adjustment cost")
 m <= parameter(:lamw, 1.5, description = "# wage markup")
 m <= parameter(:ϕh , 2., description = "# inverse frisch elasticity")
 m <= parameter(:Φw , 10., description = "# rotemberg cost for wages")
 m <= parameter(:lamf, 1.5 , description = "# price markup")
 m <= parameter(:Φp , 1., description = "# rotemberg cost for prices")
 m <= parameter(:ρR , 0.75, description = "# persistence in taylor rule")
 m <= parameter(:ψπ , 1.5, description = "# weight on inflation in taylor rule")
 m <= parameter(:ψy , 0.5, description = "# weight on output growth in taylor rule")

    # Setting steady-state parameters
    nx = get_setting(m, :nx)
    ns = get_setting(m, :ns)

    # Steady state grids for functional/distributional variables and market-clearing discount rate
    m <= SteadyStateParameterGrid(:lstar, fill(NaN, nx*ns), description = "Steady-state expected discounted
                                  marginal utility of consumption", tex_label = "l_*")
    m <= SteadyStateParameterGrid(:cstar, fill(NaN, nx*ns), description = "Steady-state consumption",
                                  tex_label = "c_*")
    m <= SteadyStateParameterGrid(:μstar, fill(NaN, nx*ns), description = "Steady-state cross-sectional
                                  density of cash on hand", tex_label = "\\mu_*")
    m <= SteadyStateParameter(:βstar, NaN, description = "Steady-state discount factor",
                              tex_label = "\\beta_*")
end

"""
```
aggregate_steadystate!(m::HetDSGE)
```
Steady state of aggregate scalar variables (Break out the "analytic" steady-state solution here instead of in steadystate!() since some of these parameters are required to construct the grids that are used in the functional/distributional steady state variable calculation (where init_grids! is called before steadystate! in model construction).
"""
function aggregate_steadystate!(m::HetDSGE)
    m <= SteadyStateParameter(:Rkstar, m[:r] + m[:δ], description = "Rental rate on capital", tex_label = "Rk_*")
    m <= SteadyStateParameter(:ωstar, (m[:α]^(m[:α]/(1-m[:α])))*(1-m[:α])*m[:Rkstar]^(-m[:α]/(1-m[:α])), description = "Real wage", tex_label = "\\omega_*")
    m <= SteadyStateParameter(:klstar, (m[:α]/(1-m[:α]))*(m[:ωstar]/m[:Rkstar])*exp(m[:γ]), description = "Capital/Labor ratio", tex_label = "kl_*")
    m <= SteadyStateParameter(:kstar, m[:klstar]*m[:H], description = "Capital", tex_label = "k_*")
    m <= SteadyStateParameter(:xstar, (1-(1-m[:δ])*exp(-m[:γ]))*m[:kstar], description = "Investment", tex_label = "x_*")
    m <= SteadyStateParameter(:ystar, (exp(-m[:α]*m[:γ])*(m[:kstar]^m[:α])*m[:H]^(1-m[:α])), description = "GDP", tex_label = "y_*")
    m <= SteadyStateParameter(:Tstar, m[:Rkstar]*m[:kstar]*exp(-m[:γ]) - m[:xstar] - (1-(1/m[:g]))*m[:ystar], description = "Net transfer to households", tex_label = "T_*")

    return m
end

"""
```
init_grids!(m::HetDSGE)
```
"""
function init_grids!(m::HetDSGE)


    nx      = get_setting(m, :nx)
    ns      = get_setting(m, :ns)
    λ       = get_setting(m, :λ)

    grids = OrderedDict()

    # Skill grid
    #lsgrid, sprob, sscale = tauchen86(m[:μ_sp].value, m[:ρ_sp].value, m[:σ_sp].value, ns, λ)
    sprob = [[0.9 0.1];[0.1 0.9]]
    sgrid = [0.8;1.2]
    (λs, vs) = eigen(Matrix{Float64}(sprob'))
    order_λs = sortperm(λs, rev = true)
    vs = vs[:,order_λs]
    ss_skill_distr = vs[:,1]/sum(vs[:,1])
    means = ss_skill_distr'*sgrid
    meanz = (get_setting(m, :zhi)+get_setting(m, :zlo))/2.
    sgrid = sgrid/(meanz*means) # so that skills integrate to 1
    sscale = sgrid[2] - sgrid[1]
    swts = (sscale/ns)*ones(ns)
    #sgrid = exp.(lsgrid)
    grids[:sgrid] = Grid(sgrid, swts, sscale)

 # Calculate endogenous grid bounds (bounds are analytic functions of model parameters,
    # which ensure there is no mass in the distributions across x where there shouldn't be)
    # Calculate the lowest possible cash on hand in steady state
    zlo = get_setting(m, :zlo)
    smin = minimum(sgrid)*zlo #exp(m[:μ_sp]/(1-m[:ρ_sp]) - get_setting(m, :λ)*sqrt(m[:σ_sp]^2/(1-m[:ρ_sp])^2))*get_setting(m, :zlo)
    m <= Setting(:xlo, m[:ωstar]*smin*m[:H] - (1+m[:r])*m[:η]*exp(-m[:γ]) + m[:Tstar] + sgrid[1]*m[:ωstar]*m[:H]*0.05, "Lower bound on cash on hand")
    m <= Setting(:xhi, max(get_setting(m, :xlo)*2, get_setting(m, :xlo)+5.), "Upper Bound on cash on hand")
    m <= Setting(:xscale, get_setting(m, :xhi) - get_setting(m, :xlo), "Size of the xgrid")
    xlo     = get_setting(m, :xlo)
    xhi     = get_setting(m, :xhi)
    xscale  = get_setting(m, :xscale)

    # Cash on hand grid
    grids[:xgrid] = Grid(uniform_quadrature(xscale), xlo, xhi, nx, scale = xscale)


    # Markov transition matrix for skill
    grids[:fgrid] = sprob./swts

    # Total grid vectorized across both dimensions
    grids[:sgrid_total] = kron(sgrid, ones(nx))
    grids[:xgrid_total] = kron(ones(ns), grids[:xgrid].points)
    grids[:weights_total] = kron(swts, grids[:xgrid].weights)

    m.grids = grids
end

function model_settings!(m::HetDSGE)
    default_settings!(m)

    # Defaults
    # Data settings for released and conditional data. Default behavior is to set vintage
    # of data to today's date.
    vint = Dates.format(now(), DSGE_DATE_FORMAT)
    m <= Setting(:data_vintage, vint, true, "vint", "Data vintage")

    saveroot = normpath(joinpath(dirname(@__FILE__), "../../../","save"))
    datapath = normpath(joinpath(dirname(@__FILE__), "../../../","save","input_data"))

    m <= Setting(:saveroot, saveroot, "Root of data directory structure")
    m <= Setting(:dataroot, datapath, "Input data directory path")

    # Solution method
    m <= Setting(:solution_method, :klein)

    m <= Setting(:krylov_reduce, false)

    # Likelihood method
    m <= Setting(:use_chand_recursion, true)

    # Anticipated shocks
    m <= Setting(:n_anticipated_shocks, 0,
                 "Number of anticipated policy shocks")
    # m <= Setting(:n_anticipated_shocks_padding, 20,
                 # "Padding for anticipated policy shocks")

    # Whether one wishes to re-compute βstar
    m <= Setting(:use_last_βstar, false, "Flag to avoid recomputing steady-state.")

    # Number of states and jumps
    m <= Setting(:normalize_distr_variables, true, "Whether or not to perform the
                 normalization of the distributional states in the Klein solution step")
    m <= Setting(:n_predetermined_variables, 0, "Number of predetermined variables after
                 removing the densities. Calculated with the Jacobian normalization.
                 This set to 0 is just the default setting, since it will always be
                 overwritten once the Jacobian is calculated.")

    m <= Setting(:state_indices, 1:3, "Which indices of m.endogenous_states correspond to
                 backward looking state variables")
    m <= Setting(:jump_indices, 4:9, "Which indices of m.endogenous_states correspond to jump
                 variables")

    # Mollifier setting parameters
    m <= Setting(:In, 0.443993816237631, "Normalizing constant for the mollifier")
    m <= Setting(:zlo, 1.0, "Lower bound on second income shock to mollify actual income")
    m <= Setting(:zhi, 5.0, "Upper bound on second income shock to mollify actual income")

    # s: Skill Distribution/ "Units of effective labor" Grid Setup
    m <= Setting(:ns, 2, "Skill distribution grid points")
    m <= Setting(:λ, 2.0, "The λ parameter in the Tauchen distribution calculation")

    # x: Cash on Hand Grid Setup
    m <= Setting(:nx, 300, "Cash on hand distribution grid points")

    m <= Setting(:trunc_distr, true)
    m <= Setting(:rescale_weights, true)
    m <= Setting(:mindens, 1e-8)

    # Function-valued variables include distributional variables
    m <= Setting(:n_function_valued_backward_looking_states, 1, "Number of function-valued
                 backward looking state variables")
    m <= Setting(:n_backward_looking_distributional_vars, 1, "Number of state variables that are
                 distributional variables.")
    m <= Setting(:n_function_valued_jumps, 1, "Number of function-valued jump variables")
    m <= Setting(:n_jump_distributional_vars, 0, "Number of jump variables that are distributional
                 variables.")

    # Note, these settings assume normalization.
    # The n degrees of freedom removed depends on the distributions/dimensions
    # of heterogeneity that we have discretized over, in this case,
    # cash on hand and the skill distribution. In general the rule of
    # thumb is, remove one degree of freedom for the first endogenous distribution (cash on
    # hand), then one additional degree of freedom for each exogenous distribution (skill
    # distribution). Multiple endogenous distributions only permit removing a single degree
    # of freedom since it is then non-trivial to obtain the marginal distributions.
    m <= Setting(:n_degrees_of_freedom_removed, 2, "Number of degrees of freedom from the
                 distributional variables to remove.")
    n_dof_removed = get_setting(m, :n_degrees_of_freedom_removed)

    ####################################################
    # Calculating the number of backward-looking states
    ####################################################
    n_backward_looking_distr_vars = get_setting(m, :n_backward_looking_distributional_vars)
    m <= Setting(:backward_looking_states_normalization_factor,
                 n_dof_removed*n_backward_looking_distr_vars, "The number of dimensions removed from the
                 backward looking state variables for the normalization.")

    n = get_setting(m, :nx)*get_setting(m, :ns)
    n_backward_looking_vars = length(get_setting(m, :state_indices))
    n_backward_looking_function_valued_vars = get_setting(m, :n_function_valued_backward_looking_states)
    n_backward_looking_scalar_vars = n_backward_looking_vars - n_backward_looking_function_valued_vars

    m <= Setting(:n_backward_looking_states, n*n_backward_looking_distr_vars +
                 n_backward_looking_scalar_vars - get_setting(m, :backward_looking_states_normalization_factor),
                 "Number of state variables, in the true sense (fully
                  backward looking) accounting for the discretization across the grid
                  of function-valued variables and the normalization of
                  distributional variables.")

    ##################################
    # Calculating the number of jumps
    ##################################
    n_jump_distr_vars = get_setting(m, :n_jump_distributional_vars)
    m <= Setting(:jumps_normalization_factor,
                 n_dof_removed*n_jump_distr_vars, "The number of dimensions removed from the
                 jump variables for the normalization.")

    n_jump_vars = length(get_setting(m, :jump_indices))
    n_jump_function_valued_vars = get_setting(m, :n_function_valued_jumps)
    n_jump_scalar_vars = n_jump_vars - n_jump_function_valued_vars

    m <= Setting(:n_jumps, n*n_jump_function_valued_vars +
                 n_jump_scalar_vars - get_setting(m, :jumps_normalization_factor),
                 "Number of jump variables (forward looking) accounting for
                  the discretization across the grid of function-valued variables and the
                  normalization of distributional variables.")

    m <= Setting(:n_model_states, get_setting(m, :n_backward_looking_states) + get_setting(m, :n_jumps),
                 "Number of 'states' in the state space model. Because backward and forward
                 looking variables need to be explicitly tracked for the Klein solution
                 method, we have n_states and n_jumps")
end

function setup_indices!(m::HetDSGE)
    nx = get_setting(m, :nx)
    ns = get_setting(m, :ns)
    endo = m.endogenous_states
    eqconds = m.equilibrium_conditions
    nxns = nx*ns

    # Endogenous function-valued states
    endo[:l′_t1]    = 1:nxns            # lagged ell function
    endo[:μ′_t1]    = nxns+1:2*nxns     # lagged wealth distribution
    #endogenous scalar-valued states
    endo[:k′_t]  = 2*nxns+1:2*nxns+1             # capital –dont get confused with steadystate object K
    endo[:R′_t1] = 2*nxns+2:2*nxns+2              # lagged real interest rate
    endo[:i′_t1] = 2*nxns+3:2*nxns+3             # lagged nominal interest rate
    endo[:π′_t1] = 2*nxns+4:2*nxns+4             # lagged inflation
    endo[:π′_t2] = 2*nxns+5:2*nxns+5             # double-lagged inflation
    endo[:π′_t3] = 2*nxns+6:2*nxns+6             # triple-lagged inflation
    endo[:y′_t1] = 2*nxns+7:2*nxns+7             # lagged gdp
    endo[:y′_t2] = 2*nxns+8:2*nxns+8             # double-lagged gdp
    endo[:y′_t3] = 2*nxns+9:2*nxns+9             # triple-lagged gdp
    endo[:y′_t4] = 2*nxns+10:2*nxns+10           # quad-lagged gdp
    endo[:z′_t1] = 2*nxns+11:2*nxns+11           # lag tfp growth
    endo[:z′_t2] = 2*nxns+12:2*nxns+12           # double-lag tfp growth
    endo[:z′_t3] = 2*nxns+13:2*nxns+13           # triple-lag tfp growth
    endo[:w′_t1] = 2*nxns+14:2*nxns+14           # lag real wages
    endo[:I′_t1] = 2*nxns+15:2*nxns+15           # lag investment–don't get this confused with i, the nominal interest rate
    # exogenous scalar-valued states:
    endo[:BP]    = 2*nxns+16:2*nxns+16        # discount factor shock
    endo[:GP]    = 2*nxns+17:2*nxns+17        # govt spending
    endo[:ZP]    = 2*nxns+18:2*nxns+18        # tfp growth
    endo[:MUP]   = 2*nxns+19:2*nxns+19        # investment shock
    endo[:LAMWP] = 2*nxns+20:2*nxns+20        # wage markup
    endo[:LAMFP] = 2*nxns+21:2*nxns+21        # price markup
    endo[:MONP]  = 2*nxns+22:2*nxns+22        # monetary policy shock

 # function-valued jumps
    endo[:l′_t]  = 2*nxns+23:3*nxns+22 # ell function
    endo[:μ′_t]  = 3*nxns+23:4*nxns+22 # wealth distribution
    #scalar-valued jumps
    endo[:R′_t]  = 4*nxns+23:4*nxns+23        # real interest rate
    endo[:i′_t]  = 4*nxns+24:4*nxns+24        # nominal interest rate
    endo[:t′_t]  = 4*nxns+25:4*nxns+25        # transfers + dividends
    endo[:w′_t]  = 4*nxns+26:4*nxns+26        # real wage
    endo[:L′_t]  = 4*nxns+27:4*nxns+27        # hours worked
    endo[:π′_t]  = 4*nxns+28:4*nxns+28        # inflation
    endo[:wageinflation′_t] = 4*nxns+29:4*nxns+29        # nominal wage inflation
    endo[:mu′_t] =  4*nxns+30:4*nxns+30        # average marginal utility
    endo[:y′_t]  = 4*nxns+31:4*nxns+31       # gdp
    endo[:I′_t]  = 4*nxns+32:4*nxns+32        # investment
    endo[:mc′_t] = 4*nxns+33:4*nxns+33        # marginal cost - this is ζ in HetDSGEₖd.pdf
    endo[:Q′_t]  = 4*nxns+34:4*nxns+34        # Tobin's qfunction
    endo[:capreturn′_t] = 4*nxns+35:4*nxns+35        # return on capital

    nvars = 4*nxns+35
    nscalars = 35 # num eqs which output scalars
    nyscalars = 13 # num scalar jumps
    nxscalars = nscalars - nyscalars # num scalar states

endo[:l_t1] = endo[:l′_t1] .+ nvars #LELL
endo[:μ_t1] = endo[:μ′_t1] .+ nvars #LM
endo[:k_t] = endo[:k′_t] .+ nvars #KK
endo[:R_t1] = endo[:R′_t1] .+ nvars #LRR
endo[:i_t1] = endo[:i′_t1] .+ nvars #LII
endo[:π_t1] = endo[:π′_t1]  .+ nvars #LPI
endo[:π_t2] = endo[:π′_t2]  .+ nvars #L2PI
endo[:π_t3] = endo[:π′_t3]  .+ nvars #L3PI
endo[:y_t1] = endo[:y′_t1] .+ nvars #LY
endo[:y_t2] = endo[:y′_t2] .+ nvars #L2Y
endo[:y_t3] = endo[:y′_t3] .+ nvars #l3Y
endo[:y_t4] = endo[:y′_t4] .+ nvars #L4Y
endo[:z_t1] = endo[:z′_t1] .+ nvars #LZ
endo[:z_t2] = endo[:z′_t2] .+ nvars #L2Z
endo[:z_t3] = endo[:z′_t3] .+ nvars #L3Z
endo[:w_t1] = endo[:w′_t1] .+ nvars #LW
endo[:I_t1] = endo[:I′_t1] .+ nvars #LX
endo[:B] = endo[:BP] .+ nvars #B
endo[:G] = endo[:GP] .+ nvars #G
endo[:Z] = endo[:ZP] .+ nvars #Z
endo[:MU] = endo[:MUP] .+ nvars #MU
endo[:LAMW] = endo[:LAMWP] .+ nvars #LAMW
endo[:LAMF] = endo[:LAMFP] .+ nvars #LAMF
endo[:MON] = endo[:MONP] .+ nvars #MON
endo[:ELL] = endo[:l′_t] .+ nvars #ELL
endo[:μ_t] = endo[:μ′_t] .+ nvars #M
endo[:R_t] = endo[:R′_t] .+ nvars #RR
endo[:i_t] = endo[:i′_t] .+ nvars #II
endo[:t_t] = endo[:t′_t]  .+ nvars #TT
endo[:w_t] = endo[:w′_t] .+ nvars #W
endo[:L_t] = endo[:L′_t] .+ nvars #HH
endo[:π_t] = endo[:π′_t] .+ nvars #PI
endo[:wageinflation_t] = endo[:wageinflation′_t] .+ nvars #PIW
endo[:mu_t] = endo[:mu′_t] .+ nvars #LAM
endo[:y_t] = endo[:y′_t] .+ nvars #Y
endo[:I_t] = endo[:I′_t]  .+ nvars #X
endo[:mc_t] = endo[:mc′_t]  .+ nvars #MC
endo[:Q_t] = endo[:Q′_t]  .+ nvars #Q
endo[:capreturn_t] = endo[:capreturn′_t] .+ nvars #RK

# create objects needed for solve.jl
# we will order function blocks as follows:
# 1. all function blocks which output a function (first real eqs, then lags)
# 2. all function blocks which map functions to scalars
# 3. all scalar blocks involving endogenous vbls (first real eqs, then lags)
# 4. shock processes
funops = 1:4 # which operators output a function


    # function blocks which output a function
    eqconds[:eq_euler]              = 1:nxns
    eqconds[:eq_kolmogorov_fwd]     = nx*ns+1:2*nx*ns
    eqconds[:eq_lag_ell]            = 2*nxns+1:3*nxns
    eqconds[:eq_lag_wealth]         = 3*nxns+1:4*nxns
    # function blocks which map functions to scalars
    eqconds[:eq_market_clearing]    = 4*nxns+1:4*nxns+1
    eqconds[:eq_lambda]             = 4*nxns+2:4*nxns+2
    #scalar blocks involving endogenous variables
    eqconds[:eq_transfers]            = 4*nxns+3:4*nxns+3 # transfers
    eqconds[:eq_investment]          = 4*nxns+4:4*nxns+4 # investment
    eqconds[:eq_tobin_q]              = 4*nxns+5:4*nxns+5 # tobin's q
    eqconds[:eq_capital_accumulation] = 4*nxns+6:4*nxns+6 # capital accumulation
    eqconds[:eq_wage_phillips] = 4*nxns+7:4*nxns+7 # wage phillips curve
    eqconds[:eq_price_phillips] = 4*nxns+8:4*nxns+8 # price phillips curve
    eqconds[:eq_marginal_cost]  = 4*nxns+9:4*nxns+9 # marginal cost
    eqconds[:eq_gdp]  = 4*nxns+10:4*nxns+10 # gdp
    eqconds[:eq_optimal_kl] = 4*nxns+11:4*nxns+11 # optimal K/L ratio
    eqconds[:eq_taylor] = 4*nxns+12:4*nxns+12 # taylor rule
    eqconds[:eq_fisher] = 4*nxns+13:4*nxns+13 # fisher eqn
    eqconds[:eq_nominal_wage_inflation] = 4*nxns+14:4*nxns+14 # nominal wage inflation
    # lagged variables
    eqconds[:LR] = 4*nxns+15:4*nxns+15 # LR
    eqconds[:LI] = 4*nxns+16:4*nxns+16 # LI
    eqconds[:LPI] = 4*nxns+17:4*nxns+17 # LPI
    eqconds[:L2PI] = 4*nxns+18:4*nxns+18 # L2PI
    eqconds[:L3PI] = 4*nxns+19:4*nxns+19 # L3PI
    eqconds[:LY] = 4*nxns+20:4*nxns+20 # LY
    eqconds[:L2Y] = 4*nxns+21:4*nxns+21 # L2Y
    eqconds[:L3Y] = 4*nxns+22:4*nxns+22 # L3Y
    eqconds[:L4Y] = 4*nxns+23:4*nxns+23 # L4Y
    eqconds[:LZ]  = 4*nxns+24:4*nxns+24 # LZ
    eqconds[:L2Z] = 4*nxns+25:4*nxns+25 # L2Z
    eqconds[:L3Z] = 4*nxns+26:4*nxns+26 # L3Z
    eqconds[:LW]  = 4*nxns+27:4*nxns+27 # LW
    eqconds[:LX] = 4*nxns+28:4*nxns+28 # LX
    # shocks
    eqconds[:F33] = 4*nxns+29:4*nxns+29 # discount factor B
    eqconds[:F34] = 4*nxns+30:4*nxns+30 # govt spending G
    eqconds[:F35] = 4*nxns+31:4*nxns+31 # tfp growth Z
    eqconds[:F36] = 4*nxns+32:4*nxns+32 # investment MU
    eqconds[:F37] = 4*nxns+33:4*nxns+33 # wage mkup LAMW
    eqconds[:F38] = 4*nxns+34:4*nxns+34 # price mkup LAMF
    eqconds[:F39] = 4*nxns+35:4*nxns+35 # monetary policy MON

# Total grid x*s
m <= Setting(:n, get_setting(m, :nx) * get_setting(m, :ns), "Total grid size, multiplying
                     across grid dimensions.")

m <= Setting(:nvars, 4*get_setting(m, :n) +35, "num variables")
m <= Setting(:nscalars, 35, " # num eqs which output scalars")
m <= Setting(:nyscalars, 13, "num scalar jumps")
m <= Setting(:nxscalars, 35 - 13, "num scalar states")

end
