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
    endogenous_states_augmented::OrderedDict{Symbol, Int}
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
    endogenous_states = collect([:kf′_t, :k′_t, :R′_t1, :i′_t1,
                                 :y′_t1, :w′_t1,:I′_t1,:b′_t,:g′_t,:z′_t,:μ′_t,:λ_w′_t,
                                 :λ_f′_t,:rm′_t,:l′_t,:μ′_t,:R′_t,:i′_t,:t′_t,:w′_t,:L′_t,:π′_t,
                                 :π_w′_t,:mu′_t,:y′_t,:I′_t,:mc′_t,:Q′_t,:capreturn′_t])

    # Exogenous shocks
    exogenous_shocks = collect([:b_sh,:g_sh,:z_sh,:μ_sh,:λ_w_sh, :λ_f_sh,:rm_sh])

    # Equilibrium conditions
    equilibrium_conditions = collect([:eq_euler,:eq_kolmogorov_fwd,
                                      :eq_market_clearing,:eq_lambda, :eq_transfers,
                                      :eq_investment,:eq_tobin_q,:eq_capital_accumulation,
                                      :eq_wage_phillips,:eq_price_phillips,:eq_marginal_cost,
                                      :eq_gdp,:eq_optimal_kl,:eq_taylor, :eq_fisher,
                                      :eq_nominal_wage_inflation,:LR,:LI,:LY,
                                      :LW,:LX,:eq_b,:eq_g,:eq_z,:eq_μ,:eq_λ_w,
                                      :eq_λ_f, :eq_rm])

    # Additional states added after solving model
    # Lagged states and observables measurement error
    endogenous_states_augmented = [:i_t1]

    # Observables
    observables = keys(m.observable_mappings)

    ########################################################################################
    # Setting indices of endogenous_states and equilibrium conditions manually for now

    setup_indices!(m)
    endo = m.endogenous_states_unnormalized
    eqcond = equilibrium_conditions

    ########################################################################################

    m.normalized_model_states = [:kf′_t]

    for (i,k) in enumerate(exogenous_shocks);            m.exogenous_shocks[k]            = i end
    for (i,k) in enumerate(endogenous_states_augmented); m.endogenous_states_augmented[k] = i+length(endogenous_states) end
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
            OrderedDict{Symbol,UnitRange}(),
            OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(),

            spec,
            subspec,
            settings,
            test_settings,
            rng,
            testing,
            OrderedDict{Symbol,Observable}())

    m <= Setting(:nx, 300)
    m <= Setting(:ns, 2)

    default_settings!(m)

    # # Set observable transformations
    init_observable_mappings!(m)

    # Initialize model indices
    init_model_indices!(m)

    # Set settings
    model_settings!(m)
    # default_test_settings!(m)
    for custom_setting in values(custom_settings)
        m <= custom_setting
    end

    # Initialize parameters
    init_parameters!(m)

    # Initialize aggregate steady state parameters (necessary for grid construction)
    aggregate_steadystate!(m)

    # Initialize grids
    init_grids!(m)


    # # Solve for the steady state
    steadystate!(m)

    # # So that the indices of m.endogenous_states reflect the normalization
    normalize_model_state_indices!(m)

    init_subspec!(m)

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
    ######################################
    # Parameters that affect steady-state
    ######################################
    m <= parameter(:r, 0.01, fixed = true,
                   description = "r: Steady-state real interest rate.", tex_label = "r")
    m <= parameter(:α, 0.3, fixed = true, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), Normal(0.30, 0.05),
                   description = "α: Capital elasticity in the intermediate goods
                   sector's production function (also known as the capital share).",
                   tex_label = "\\alpha")
        # Check this: Previously the description was "Aggregate hours worked"
        m <= parameter(:H, 1.0, (-1000., 1000.), (-1e3, 1e3), Untransformed(),
                       Normal(-45., 5.), fixed = true,
                       description = "Lmean: Mean level of hours", tex_label = "\\bar{L}")
    m <= parameter(:δ, 0.03, fixed = true,
                   description = "δ: The capital depreciation rate", tex_label = "\\delta")

    m <= parameter(:μ_sp, 0.0, fixed = true, description = "μ_sp: The trend in the skill process",
                   tex_label = "\\mu_{sp}")
    m <= parameter(:ρ_sp, 0.7, (1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = true,
                   description="ρ_sp: AR(1) coefficient in the skill process.",
                   tex_label="\\rho_{sp}")
    m <= parameter(:σ_sp, 0.01, (1e-8, 5.), (1e-8, 5.), Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description="σ_sp: The standard deviation of the skill process.",
                   tex_label="\\sigma_{sp}")

    # Exogenous processes - level
        # Uncomment scaling once adjusted properly in the code
        m <= parameter(:γ, 0.0, (-5.0, 5.0), (-5., 5.), Untransformed(),
                       Normal(0.4, 0.1), fixed = true, # scaling = x -> x/100,
                       description = "γ: The log of the steady-state growth rate of technology",
                       tex_label="100\\gamma")
    # Check this
    m <= parameter(:g, 1/(1-0.01), fixed = true,
                   description = "g_star: 1 - (c_star + i_star)/y_star",
                   tex_label = "g_*")
    # Not in m1002
    m <= parameter(:η, 0.1, description = "η: Borrowing constraint (normalized by TFP)")

    ####################################################
    # Parameters that affect dynamics (not steady-state)
    ####################################################
    m <= parameter(:spp, 4., (-15., 15.), (-15., 15.), Untransformed(),
                   Normal(4., 1.5), fixed = false,
                   description = "S'': The second derivative of households' cost of adjusting investment.",
                   tex_label = "S''")

    m <= parameter(:ϕh, 2., (1e-5, 10.), (1e-5, 10.), Exponential(),
                   Normal(2, 0.75), fixed = false, description = "inverse frisch elasticity",
                   tex_label = "\\phi_h")

    m <= parameter(:κ_p, 0.5, (1e-5, 5.), (1e-5, 5.), SquareRoot(), GammaAlt(0.5, 0.3), fixed = false, description = "κ_p : The slope of the Price Phillips curve", tex_label = "\\kappa")
    m <= parameter(:κ_w, 0.5, (1e-5, 5.), (1e-5, 5.), SquareRoot(), GammaAlt(0.5, 0.3), fixed = false, description = "κ_w: The slope of the Wage Phillips curve", tex_label = "\\kappa")

    m <= parameter(:ρR , 0.75, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(),
                   BetaAlt(0.75, 0.10), fixed = false,
                   description = "ρ: The degree of inertia in the monetary policy rule.",
                   tex_label="\\rho_R")
    m <= parameter(:ψπ , 1.5, (1e-5, 10.), (1e-5, 10.0), Exponential(),
                   Normal(1.5, 0.25), fixed = false,
                   description = "ψ1: Weight on inflation gap in monetary policy rule.",
                   tex_label = "\\psi_1")
    m <= parameter(:ψy , 0.5, (-0.5, 0.5), (-0.5, 0.5), Untransformed(),
                   Normal(0.12, 0.05), fixed = false,
                   description = "ψy: Weight on output growth in monetary policy rule")


    # Exogenous processes - autocorrelation
    m <= parameter(:ρ_G, 0.5, (1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = false,
                   description = "ρ_g: AR(1) coefficient in the government spending process.",
                   tex_label = "\\rho_g")
    m <= parameter(:ρ_B, 0.5, (1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = false,
                   description = "ρ_b: AR(1) coefficient in the intertemporal preference shifter process.",
                   tex_label = "\\rho_B")
    m <= parameter(:ρ_μ, 0.5, (1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = false,
                   description = "ρ_μ: AR(1) coefficient in capital adjustment cost process.",
                   tex_label = "\\rho_{\\mu}")
    m <= parameter(:ρ_z, 0.5,(1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = false,
                   description = "ρ_z: AR(1) coefficient in the technology process.",
                   tex_label = "\\rho_z")
    m <= parameter(:ρ_lamf, 0.5, (1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = false,
                   description = "ρ_λ_f: AR(1) coefficient in the price mark-up shock process.",
                   tex_label = "\\rho_{\\lambda_f}")
    m <= parameter(:ρ_lamw, 0.5, (1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = false,
                   description = "ρ_λ_w: AR(1) coefficient in the wage mark-up shock process.",
                   tex_label = "\\rho_{\\lambda_w}")
    m <= parameter(:ρ_mon, 0.5, (1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = false,
                   description = "ρ_rm: AR(1) coefficient in the monetary policy shock process.",
                   tex_label = "\\rho_{r^m}")

    m <= parameter(:σ_g, 0.15, (1e-8, 5.), (1e-8, 5.), Exponential(),
                   RootInverseGamma(2, 0.10), fixed = false,
                   description = "σ_g: The standard deviation of the government spending process.",
                   tex_label = "\\sigma_{g}")
    m <= parameter(:σ_b, 0.15, (1e-8, 5.), (1e-8, 5.), Exponential(),
                   RootInverseGamma(2, 0.10), fixed = false,
                   description = "σ_b: The standard deviation of the intertemporal preference shifter process.",
                   tex_label = "\\sigma_{b}")
    m <= parameter(:σ_μ, 0.15, (1e-8, 5.), (1e-8, 5.), Exponential(),
                   RootInverseGamma(2, 0.10), fixed = false,
                   description = "σ_μ: The standard deviation of the exogenous marginal efficiency of investment shock process.",
                   tex_label = "\\sigma_{\\mu}")
    m <= parameter(:σ_z, 0.15,(1e-8, 5.), (1e-8, 5.), Exponential(),
                   RootInverseGamma(2, 0.10), fixed = false,
                   description = "σ_z: The standard deviation of the process describing the stationary component of productivity.",
                   tex_label = "\\sigma_z")
    m <= parameter(:σ_λ_f, 0.15, (1e-8, 5.), (1e-8, 5.), Exponential(),
                   RootInverseGamma(2, 0.10), fixed = false,
                   description = "σ_λ_f: The mean of the process that generates the price elasticity of the composite good. Specifically, the elasticity is (1+λ_{f,t})/(λ_{f,t}).",
                   tex_label = "\\sigma_{\\lambda_f}")
    m <= parameter(:σ_λ_w, 0.15, (1e-8, 5.), (1e-8, 5.), Exponential(),
                   RootInverseGamma(2, 0.10), fixed = false,
                   description = "σ_λ_w",
                   tex_label = "\\sigma_{\\lambda_w}")
    m <= parameter(:σ_rm, 0.15, (1e-8, 5.), (1e-8, 5.), Exponential(),
                   RootInverseGamma(2, 0.10), fixed = false,
                   description = "σ_r_m: The standard deviation of the monetary policy shock.", tex_label = "\\sigma_{r^m}")

    m <= parameter(:π_star, 0.7000, (1e-5, 10.), (1e-5, 10.), Exponential(), GammaAlt(0.62, 0.1), fixed=false, scaling = x -> 1 + x/100,
               description="π_star: The steady-state rate of inflation.",
               tex_label="\\pi_*")


    m <= parameter(:Lmean, -45.9364, (-1000., 1000.), (-1e3, 1e3), Untransformed(), Normal(-45., 5.), fixed=false,
                   description="Lmean: Mean level of hours.",
                   tex_label="\\bar{L}")

m <= parameter(:e_y, 0.0, fixed = true, description = "e_y: Measurement error on GDP", tex_label = "e_y")
m <= parameter(:e_L, 0.0, fixed = true, description = "e_L: Measurement error on hours worked", tex_label = "e_L")
m <= parameter(:e_w, 0.0, fixed = true, description = "e_w: Measurement error on wages", tex_label ="e_w")
m <= parameter(:e_π, 0.0, fixed = true, description = "e_π: Measurement error on GDP deflator", tex_label = "e_π")
m <= parameter(:e_R, 0.0, fixed = true, description = "e_R: Measurement error on nominal rate of interest", tex_label = "e_R")
m <= parameter(:e_c, 0.0, fixed = true, description = "e_c: Measurement error on consumption", tex_label = "e_c")
m <= parameter(:e_i, 0.0, fixed = true, description = "e_i: Measurement error on investment", tex_label = "e_i")


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

    endo = m.endogenous_states_unnormalized
    state_indices = [endo[:kf′_t];  endo[:k′_t]; endo[:R′_t1];endo[:i′_t1];
                     endo[:y′_t1];
                     endo[:w′_t1];endo[:I′_t1];endo[:b′_t];
                     endo[:g′_t];endo[:z′_t];endo[:μ′_t];endo[:λ_w′_t]; endo[:λ_f′_t];endo[:rm′_t]]

    jump_indices = [endo[:l′_t];endo[:R′_t];endo[:i′_t];endo[:t′_t];endo[:w′_t];
                    endo[:L′_t];endo[:π′_t];endo[:π_w′_t];endo[:mu′_t];endo[:y′_t];
                    endo[:I′_t];endo[:mc′_t]; endo[:Q′_t];endo[:capreturn′_t]]

    m <= Setting(:state_indices, state_indices, "Which indices of m.endogenous_states correspond to
                 backward looking state variables")
    m <= Setting(:jump_indices, jump_indices, "Which indices of m.endogenous_states correspond to jump
                 variables")
    #m.state_variables = m.endogenous_states.keys[get_setting(m, :state_indices)]
    #m.jump_variables = m.endogenous_states.keys[get_setting(m, :jump_indices)]


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
    m <= Setting(:n_jump_distributional_vars, 1, "Number of jump variables that are distributional
                 variables.")

    # Note, these settings assume normalization.
    # The n degrees of freedom removed depends on the distributions/dimensions
    # of heterogeneity that we have discretized over, in this case,
    # cash on hand and the skill distribution. In general the rule of
    # thumb is, remove one degree of freedom for the first endogenous distribution (cash on
    # hand), then one additional degree of freedom for each exogenous distribution (skill
    # distribution). Multiple endogenous distributions only permit removing a single degree
    # of freedom since it is then non-trivial to obtain the marginal distributions.
    m <= Setting(:n_degrees_of_freedom_removed_state, 2, "Number of degrees of freedom from the distributional variables to remove.")
    m <= Setting(:n_degrees_of_freedom_removed_jump, 0, "Number of degrees of freedom from the distributional variables to remove.")
    n_dof_removed_state = get_setting(m, :n_degrees_of_freedom_removed_state)
    n_dof_removed_jump = get_setting(m, :n_degrees_of_freedom_removed_jump)

    ####################################################
    # Calculating the number of backward-looking states
    ####################################################
    n_backward_looking_distr_vars = get_setting(m, :n_backward_looking_distributional_vars)
    m <= Setting(:backward_looking_states_normalization_factor,
                 n_dof_removed_state*n_backward_looking_distr_vars, "The number of dimensions removed from the
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
                 n_dof_removed_jump*n_jump_distr_vars, "The number of dimensions removed from the
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
    m <= Setting(:n_model_states_augmented, get_setting(m, :n_model_states) +
                 length(m.endogenous_states_augmented))
end

function setup_indices!(m::HetDSGE)
    nx = get_setting(m, :nx)
    ns = get_setting(m, :ns)
    endo = m.endogenous_states_unnormalized
    eqconds = m.equilibrium_conditions
    nxns = nx*ns

    # Endogenous function-valued states
    endo[:kf′_t]    = 1:nxns    #  combination of lagged ell function and lagged m function that predicts m

    #endogenous scalar-valued states
    endo[:k′_t]  = nxns+1:nxns+1             # capital –dont get confused with steadystate object K
    endo[:R′_t1] = nxns+2:nxns+2              # lagged real interest rate
    endo[:i′_t1] = nxns+3:nxns+3             # lagged nominal interest rate
    endo[:y′_t1] = nxns+4:nxns+4             # lagged gdp
    endo[:w′_t1] = nxns+5:nxns+5           # lag real wages
    endo[:I′_t1] = nxns+6:nxns+6           # lag investment–don't get this confused with i, the nominal interest rate
    # exogenous scalar-valued states:
    endo[:b′_t]    = nxns+7:nxns+7        # discount factor shock
    endo[:g′_t]    = nxns+8:nxns+8        # govt spending
    endo[:z′_t]  = nxns+9:nxns+9        # tfp growth
    endo[:μ′_t]   = nxns+10:nxns+10        # investment shock
    endo[:λ_w′_t] = nxns+11:nxns+11        # wage markup
    endo[:λ_f′_t] = nxns+12:nxns+12        # price markup
    endo[:rm′_t]  = nxns+13:nxns+13        # monetary policy shock

    # function-valued jumps
    endo[:l′_t]  = nxns+14:2*nxns+13 # ell function

    #scalar-valued jumps
    endo[:R′_t]  = 2*nxns+14:2*nxns+14        # real interest rate
    endo[:i′_t]  = 2*nxns+15:2*nxns+15        # nominal interest rate
    endo[:t′_t]  = 2*nxns+16:2*nxns+16        # transfers + dividends
    endo[:w′_t]  = 2*nxns+17:2*nxns+17        # real wage
    endo[:L′_t]  = 2*nxns+18:2*nxns+18        # hours worked
    endo[:π′_t]  = 2*nxns+19:2*nxns+19        # inflation
    endo[:π_w′_t] = 2*nxns+20:2*nxns+20        # nominal wage inflation
    endo[:mu′_t] =  2*nxns+21:2*nxns+21        # average marginal utility
    endo[:y′_t]  = 2*nxns+22:2*nxns+22       # gdp
    endo[:I′_t]  = 2*nxns+23:2*nxns+23        # investment
    endo[:mc′_t] = 2*nxns+24:2*nxns+24        # marginal cost - this is ζ in HetDSGEₖd.pdf
    endo[:Q′_t]  = 2*nxns+25:2*nxns+25        # Tobin's qfunction
    endo[:capreturn′_t] = 2*nxns+26:2*nxns+26        # return on capital

    nvars = 2*nxns+26
    nscalars = 26 # num eqs which output scalars
    nyscalars = 13 # num scalar jumps
    nxscalars = nscalars - nyscalars # num scalar states

# create objects needed for solve.jl
# we will order function blocks as follows:
# 1. all function blocks which output a function (first real eqs, then lags)
# 2. all function blocks which map functions to scalars
# 3. all scalar blocks involving endogenous vbls (first real eqs, then lags)
# 4. shock processes
funops = 1:2 # which operators output a function


    # function blocks which output a function
    eqconds[:eq_euler]              = 1:nxns
    eqconds[:eq_kolmogorov_fwd]     = nx*ns+1:2*nx*ns

    # function blocks which map functions to scalars
    eqconds[:eq_market_clearing]    = 2*nxns+1:2*nxns+1
    eqconds[:eq_lambda]             = 2*nxns+2:2*nxns+2
    #scalar blocks involving endogenous variables
    eqconds[:eq_transfers]            = 2*nxns+3:2*nxns+3 # transfers
    eqconds[:eq_investment]          = 2*nxns+4:2*nxns+4 # investment
    eqconds[:eq_tobin_q]              = 2*nxns+5:2*nxns+5 # tobin's q
    eqconds[:eq_capital_accumulation] = 2*nxns+6:2*nxns+6 # capital accumulation
    eqconds[:eq_wage_phillips] = 2*nxns+7:2*nxns+7 # wage phillips curve
    eqconds[:eq_price_phillips] = 2*nxns+8:2*nxns+8 # price phillips curve
    eqconds[:eq_marginal_cost]  = 2*nxns+9:2*nxns+9 # marginal cost
    eqconds[:eq_gdp]  = 2*nxns+10:2*nxns+10 # gdp
    eqconds[:eq_optimal_kl] = 2*nxns+11:2*nxns+11 # optimal K/L ratio
    eqconds[:eq_taylor] = 2*nxns+12:2*nxns+12 # taylor rule
    eqconds[:eq_fisher] = 2*nxns+13:2*nxns+13 # fisher eqn
    eqconds[:eq_nominal_wage_inflation] = 2*nxns+14:2*nxns+14 # nominal wage inflation
    # lagged variables
    eqconds[:LR] = 2*nxns+15:2*nxns+15 # LR
    eqconds[:LI] = 2*nxns+16:2*nxns+16 # LI
    eqconds[:LY] = 2*nxns+17:2*nxns+17 # LY
    eqconds[:LW]  = 2*nxns+18:2*nxns+18 # LW
    eqconds[:LX] = 2*nxns+19:2*nxns+19 # LX
    # shocks
    eqconds[:eq_b] = 2*nxns+20:2*nxns+20 # discount factor B
    eqconds[:eq_g] = 2*nxns+21:2*nxns+21 # govt spending G
    eqconds[:eq_z] = 2*nxns+22:2*nxns+22 # tfp growth Z
    eqconds[:eq_μ] = 2*nxns+23:2*nxns+23 # investment MU
    eqconds[:eq_λ_w] = 2*nxns+24:2*nxns+24 # wage mkup LAMW
    eqconds[:eq_λ_f] = 2*nxns+25:2*nxns+25 # price mkup LAMF
    eqconds[:eq_rm] = 2*nxns+26:2*nxns+26 # monetary policy MON

    # Total grid x*s
    m <= Setting(:n, get_setting(m, :nx) * get_setting(m, :ns), "Total grid size, multiplying
                         across grid dimensions.")

m <= Setting(:nvars, 2*get_setting(m, :n) + 26, "num variables")
m <= Setting(:nscalars, 26, " # num eqs which output scalars")
m <= Setting(:nyscalars, 13, "num scalar jumps")
m <= Setting(:nxscalars, 26 - 13, "num scalar states")

m.endogenous_states = deepcopy(endo)

 #= endo = m.endogenous_states_unnormalized
    state_indices = [endo[:kf′_t];  endo[:k′_t]; endo[:R′_t1];endo[:i′_t1];
                     endo[:y′_t1];
                     endo[:w′_t1];endo[:I′_t1];endo[:B′];
                     endo[:G′];endo[:z′_t];endo[:MU′];endo[:LAMW′]; endo[:LAMF′];endo[:MON′]]

    jump_indices = [endo[:l′_t];endo[:R′_t];endo[:i′_t];endo[:t′_t];endo[:w′_t];
                    endo[:L′_t];endo[:π′_t];endo[:π_w′_t];endo[:mu′_t];endo[:y′_t];
                    endo[:I′_t];endo[:mc′_t]; endo[:Q′_t];endo[:capreturn′_t]]

    m <= Setting(:state_indices, state_indices, "Which indices of m.endogenous_states correspond to
                 backward looking state variables")
    m <= Setting(:jump_indices, jump_indices, "Which indices of m.endogenous_states correspond to jump
                 variables")=#


end
