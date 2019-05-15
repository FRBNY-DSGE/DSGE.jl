import DataStructures: OrderedDict
"""
```
HetDSGEGovDebt{T} <: AbstractHeterogeneousModel{T}
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
mutable struct HetDSGEGovDebt{T} <: AbstractHetModel{T}
    parameters::ParameterVector{T}               # vector of all time-invariant model parameters
    steady_state::ParameterVector{T}             # model steady-state values

    # Temporary to get it to work. Need to
    # figure out a more flexible way to define
    # "grids" that are not necessarily quadrature
    # grids within the model
    grids::OrderedDict{Symbol,Union{Grid, Array}}
    keys::OrderedDict{Symbol,Int}                    # Human-readable names for all the model
                                                     # parameters and steady-states

    state_variables::Vector{Symbol}                  # Vector of symbols of the state variables
    jump_variables::Vector{Symbol}                   # Vector of symbols of the jump variables
    normalized_model_states::Vector{Symbol}          # All of the distributional model
                                                     # state variables that need to be normalized
    # Vector of unnormalized ranges of indices
    endogenous_states_unnormalized::OrderedDict{Symbol,UnitRange}
    # Vector of ranges corresponding to normalized (post Klein solution) indices
    endogenous_states::OrderedDict{Symbol,UnitRange}
    endogenous_states_original::OrderedDict{Symbol,UnitRange}

    exogenous_shocks::OrderedDict{Symbol,Int}
    expected_shocks::OrderedDict{Symbol,Int}
    equilibrium_conditions::OrderedDict{Symbol,UnitRange}
    endogenous_states_augmented::OrderedDict{Symbol, Int}
    observables::OrderedDict{Symbol,Int}

    spec::String                                     # Model specification number (eg "m990")
    subspec::String                                  # Model subspecification (eg "ss0")
    settings::Dict{Symbol,Setting}                   # Settings/flags for computation
    test_settings::Dict{Symbol,Setting}              # Settings/flags for testing mode
    rng::MersenneTwister                             # Random number generator
    testing::Bool                                    # Whether we are in testing mode or not
    observable_mappings::OrderedDict{Symbol, Observable}
end

description(m::HetDSGEGovDebt) = "HetDSGEGovDebt, $(m.subspec)"

"""
`init_model_indices!(m::HetDSGEGovDebt)`

Arguments:
`m:: HetDSGEGovDebt`: a model object

Description:
Initializes indices for all of `m`'s states, shocks, and equilibrium conditions.
"""
function init_model_indices!(m::HetDSGEGovDebt, states::Vector{Symbol}, jumps::Vector{Symbol})

    endogenous_states = collect(vcat(states, jumps))

    # Exogenous shocks
    exogenous_shocks = collect([:b_sh,:g_sh,:z_sh,:μ_sh,:λ_w_sh, :λ_f_sh,:rm_sh])

    # Equilibrium conditions
    equilibrium_conditions = collect([:eq_euler,:eq_kolmogorov_fwd,
                                      :eq_market_clearing,:eq_lambda, :eq_transfers,
                                      :eq_investment,:eq_tobin_q,:eq_capital_accumulation,
                                      :eq_wage_phillips,:eq_price_phillips,:eq_marginal_cost,
                                      :eq_gdp,:eq_optimal_kl,:eq_taylor, :eq_fisher,
                                      :eq_nominal_wage_inflation, :eq_fiscal_rule,
                                      :eq_g_budget_constraint, :LR,:LI,:LY,
                                      :LW,:LX,:eq_b,:eq_g,:eq_z,:eq_μ,:eq_λ_w,
                                      :eq_λ_f, :eq_rm]) #, :eq_consumption])

    # Additional states added after solving model
    # Lagged states and observables measurement error
    endogenous_states_augmented = [:i_t1, :c_t, :c_t1]

    # Observables
    observables = keys(m.observable_mappings)

    ########################################################################################
    # Setting indices of endogenous_states and equilibrium conditions manually for now

    setup_indices!(m)
    endo = m.endogenous_states_unnormalized
    m.endogenous_states_original = deepcopy(endo)
    eqcond = equilibrium_conditions
    ########################################################################################

    m.normalized_model_states = [:kf′_t]

    for (i,k) in enumerate(exogenous_shocks); m.exogenous_shocks[k] = i end
    for (i,k) in enumerate(observables);      m.observables[k]      = i end
end

function HetDSGEGovDebt(subspec::String="ss0";
                   custom_settings::Dict{Symbol, Setting} = Dict{Symbol, Setting}(),
                   testing = false, testing_gamma::Bool = false)

    # Model-specific specifications
    spec               = "het_dsge"
    subspec            = subspec
    settings           = Dict{Symbol,Setting}()
    test_settings      = Dict{Symbol,Setting}()
    rng                = MersenneTwister(0)

    # initialize empty model
    m = HetDSGEGovDebt{Float64}(
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
            OrderedDict{Symbol,UnitRange}(), OrderedDict{Symbol,UnitRange}(),
            OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(),

            spec,
            subspec,
            settings,
            test_settings,
            rng,
            testing,
            OrderedDict{Symbol,Observable}())

    default_settings!(m)

    # # Set observable transformations
    init_observable_mappings!(m)

    # Set settings
    model_settings!(m)
    # default_test_settings!(m)
    for custom_setting in values(custom_settings)
        m <= custom_setting
    end

    # Endogenous states
    states = collect([:kf′_t,:k′_t, :R′_t1,:i′_t1, :y′_t1,:w′_t1,:I′_t1, :bg′_t,
                      :b′_t,:g′_t,:z′_t,:μ′_t,:λ_w′_t, :λ_f′_t,:rm′_t])

    jumps = collect([:l′_t,:R′_t,:i′_t,:t′_t,:w′_t, :L′_t,:π′_t,:π_w′_t,:margutil′_t,:y′_t, :I′_t,
                     :mc′_t,:Q′_t,:capreturn′_t, :tg′_t])

    # Initialize model indices
    init_model_indices!(m, states, jumps)

    # Init what will keep track of # of states, jumps, and states/jump indices
    # (need model indices first)
    init_states_and_jumps!(m, states, jumps)

    # Initialize parameters
    init_parameters!(m, testing_gamma)

    # Initialize aggregate steady state parameters (necessary for grid construction)
    aggregate_steadystate!(m)

    # Initialize grids
    init_grids!(m)

    # Solve for the steady state
    #steadystate!(m)

    # So that the indices of m.endogenous_states reflect the normalization
    normalize_model_state_indices!(m)

    endogenous_states_augmented = [:i_t1, :c_t, :c_t1]
    for (i,k) in enumerate(endogenous_states_augmented)
        m.endogenous_states_augmented[k] = i + first(collect(values(m.endogenous_states))[end])
    end
    m <= Setting(:n_model_states_augmented, get_setting(m, :n_model_states) +
                 length(m.endogenous_states_augmented))

    init_subspec!(m)

    #=Qx, Qy, Qleft, Qright = compose_normalization_matrices(m)
    m <= Setting(:n_predetermined_variables, size(Qx, 1))
    m <= Setting(:Qleft, Qleft)
    m <= Setting(:Qright, Qright)=#
    return m
end

"""
```
init_parameters!(m::HetDSGEGovDebt)
```

Initializes the model's parameters, as well as empty values for the steady-state
parameters (in preparation for `steadystate!(m)` being called to initialize
those).
"""
function init_parameters!(m::HetDSGEGovDebt, testing_gamma::Bool)
    ######################################
    # Parameters that affect steady-state
    ######################################
    m <= parameter(:α, 0.3, fixed = false, (1e-5, 0.999), (1e-5, 0.999),
                   SquareRoot(), Normal(0.30, 0.05),
                   description = "α: Capital elasticity in the intermediate goods" *
                   "sector's production function (also known as the capital share).",
                   tex_label = "\\alpha")
    # Check this: Previously the description was "Aggregate hours worked"
    m <= parameter(:H, 1.0, (-1000., 1000.), (-1e3, 1e3), Untransformed(),
                   Normal(-45., 5.), fixed = true,
                   description = "Lmean: Mean level of hours", tex_label = "H")
    m <= parameter(:δ, 0.03, fixed = true,
                   description = "δ: The capital depreciation rate", tex_label = "\\delta")
    m <= parameter(:μ_sp, 0.0, fixed = true,
                   description = "μ_sp: The trend in the skill process",
                   tex_label = "\\mu_{sp}")

    # Exogenous processes - level
    # Uncomment scaling once adjusted properly in the code
    if testing_gamma == true
        m <= parameter(:γ, 0.0, (-5.0, 5.0), (-5., 5.), Untransformed(),
                       Normal(0.4, 0.1), fixed = false, scaling = x -> x/100,
                       description = "γ: The log of the steady-state growth rate of technology",
                       tex_label="100\\gamma")
    else
        m <= parameter(:γ, 0.5, (-5.0, 5.0), (-5., 5.), Untransformed(),
                       Normal(0.4, 0.1), fixed = false, scaling = x -> x/100,
                       description = "γ: The log of the steady-state growth rate of technology",
                       tex_label="100\\gamma")
    end

    m <= parameter(:r, 0.5, (1e-5, 10.0), (1e-5, 10.0), Exponential(),
                   GammaAlt(0.5, 0.5), fixed = false, scaling = x -> x/100,
                   description= "r: Quarterly steady-state real interest rate.",
                   tex_label= "100*r^{HetDSGE}")
    #=
    m <= parameter(:r, 0.6, (1e-5, 10.0), (1e-5, 10.), Exponential(),
                   GammaAlt(0.25, .1), fixed = false, scaling = x -> x/100 + .4/100,
                   description= "r: Quarterly steady-state real interest rate.",
                   tex_label= "100*(r^{HetDSGE}-\\gamma^{FRBNY})")
    =#

    m <= parameter(:g, 1/(1-0.01), fixed = true,
                   description = "g_star: 1 - (c_star + i_star)/y_star",
                   tex_label = "g_*")

    m <= parameter(:β_save, 0.0, fixed = true,
                   description = "saving the betas per particle",
                   tex_label = "\\beta_save")
    m <= parameter(:sH_over_sL, 6.33333, fixed = true,
                   description = "Ratio of high to low earners", tex_label = "s_H / s_L")

    m <= parameter(:pLH, 0.005, (0.0025, 0.095), (0.0025, 0.095), Untransformed(),
                   Uniform(0.005, 0.095), fixed = false,
                   description = "Prob of going from low to high persistent skill",
                   tex_label = "p(s_L \\mid s_H)")
    m <= parameter(:pHL, 0.03, (0.0025, 0.095), (0.0025, 0.095), Untransformed(),
                   Uniform(0.0025, 0.095), fixed = false,
                   description = "Prob of going from high to low persistent skill",
                   tex_label = "p(s_H \\mid s_L)")

    m <= parameter(:BoverY, 0.26, fixed = true, description = "???", tex_label = "B / Y")

    m <= parameter(:zlo, 0.0323232, fixed = true,
                   description = "Lower bound on second income shock to mollify actual income",
                   tex_label = "\\underbar{z}")
    m <= parameter(:zhi, 2-m[:zlo].value, fixed = true,
                   description = "Upper bound on second income shock to mollify actual income",
                   tex_label = "\\bar{z}")

    m <= parameter(:mpc, 0.23395, fixed = true, tex_label = "MPC")
    m <= parameter(:pc0, 0.071893, fixed = true, tex_label = "pc0")

    # Not in m1002
    m <= parameter(:η, 0.0, description = "η: Borrowing constraint (normalized by TFP)",
                   tex_label = "\\eta")

    ####################################################
    # Parameters that affect dynamics (not steady-state)
    ####################################################
    m <= parameter(:spp, 4., (-15., 15.), (-15., 15.), Untransformed(),
                   Normal(4., 1.5), fixed = false,
                   description = "S'': The second derivative of households' cost " *
                   "of adjusting investment.", tex_label = "S''")

    m <= parameter(:ϕh, 2., (1e-5, 10.), (1e-5, 10.), Exponential(),
                   Normal(2, 0.75), fixed = false, description = "inverse frisch elasticity",
                   tex_label = "\\phi_h")

    m <= parameter(:κ_p, 0.5, (1e-5, 5.), (1e-5, 5.), SquareRoot(), GammaAlt(0.5, 0.3),
                   fixed = false, description = "κ_p : The slope of the Price Phillips curve",
                   tex_label = "\\kappa_p")
    m <= parameter(:κ_w, 0.5, (1e-5, 5.), (1e-5, 5.), SquareRoot(), GammaAlt(0.5, 0.3),
                   fixed = false, description = "κ_w: The slope of the Wage Phillips curve",
                   tex_label = "\\kappa_w")

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
                   description = "ψy: Weight on output growth in monetary policy rule",
                   tex_label = "\\psi_y")

    m <= parameter(:δb, 1., (0.0, 1.0), (0.0, 1.0), Untransformed(),
                   Uniform(0.0, 1.0), fixed = false, description = "=1 means balanced budget",
                   tex_label = "\\delta_b")

    # Exogenous processes - autocorrelation
    m <= parameter(:ρ_G, 0.5, (1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = false,
                   description = "ρ_g: AR(1) coefficient in the government spending process.",
                   tex_label = "\\rho_g")
    m <= parameter(:ρ_B, 0.5, (1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = false,
                   description = "ρ_b: AR(1) coefficient in intertemporal preference " *
                   "shift process.", tex_label = "\\rho_B")
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
                   description = "σ_g: standard dev. of the government spending process.",
                   tex_label = "\\sigma_{g}")
    m <= parameter(:σ_b, 0.15, (1e-8, 5.), (1e-8, 5.), Exponential(),
                   RootInverseGamma(2, 0.10), fixed = false,
                   description = "σ_b: standard dev. of the intertemporal preference " *
                   "shifter process.", tex_label = "\\sigma_{b}")
    m <= parameter(:σ_μ, 0.15, (1e-8, 5.), (1e-8, 5.), Exponential(),
                   RootInverseGamma(2, 0.10), fixed = false,
                   description = "σ_μ: standard dev. of the exogenous marginal efficiency" *
                   " of investment shock process.", tex_label = "\\sigma_{\\mu}")
    m <= parameter(:σ_z, 0.15,(1e-8, 5.), (1e-8, 5.), Exponential(),
                   RootInverseGamma(2, 0.10), fixed = false,
                   description = "σ_z: standard dev. of the process describing the " *
                   "stationary component of productivity.", tex_label = "\\sigma_z")
    m <= parameter(:σ_λ_f, 0.15, (1e-8, 5.), (1e-8, 5.), Exponential(),
                   RootInverseGamma(2, 0.10), fixed = false,
                   description = "σ_λ_f: mean of the process that generates the price " *
                   "elasticity of the composite good. Specifically, the elasticity is " *
                   "(1+λ_{f,t})/(λ_{f,t}).", tex_label = "\\sigma_{\\lambda_f}")
    m <= parameter(:σ_λ_w, 0.15, (1e-8, 5.), (1e-8, 5.), Exponential(),
                   RootInverseGamma(2, 0.10), fixed = false,
                   description = "σ_λ_w", tex_label = "\\sigma_{\\lambda_w}")
    m <= parameter(:σ_rm, 0.15, (1e-8, 5.), (1e-8, 5.), Exponential(),
                   RootInverseGamma(2, 0.10), fixed = false,
                   description = "σ_r_m: standard dev. of the monetary policy shock.",
                   tex_label = "\\sigma_{r^m}")

    m <= parameter(:π_star, 0.7000, (1e-5, 10.), (1e-5, 10.), Exponential(),
                   GammaAlt(0.62, 0.1), fixed=false, scaling = x -> 1 + x/100,
                   description="π_star: steady-state rate of inflation.",
                   tex_label="\\pi_*")


    m <= parameter(:Lmean, -45.9364, (-1000., 1000.), (-1e3, 1e3), Untransformed(),
                   Normal(-45., 5.), fixed=false,
                   description="Lmean: Mean level of hours.", tex_label="\\bar{L}")

    m <= parameter(:e_y, 0.0, fixed = true,
                   description = "e_y: Measurement error on GDP", tex_label = "e_y")
    m <= parameter(:e_L, 0.0, fixed = true,
                   description = "e_L: Measurement error on hours worked", tex_label = "e_L")
    m <= parameter(:e_w, 0.0, fixed = true, description = "e_w: Measurement error on wages",
                   tex_label ="e_w")
    m <= parameter(:e_π, 0.0, fixed = true,
                   description = "e_π: Measurement error on GDP deflator", tex_label = "e_π")
    m <= parameter(:e_R, 0.0, fixed = true,
                   description = "e_R: Measurement error on nominal rate of interest",
                   tex_label = "e_R")
    m <= parameter(:e_c, 0.0, fixed = true,
                   description = "e_c: Measurement error on consumption", tex_label = "e_c")
    m <= parameter(:e_i, 0.0, fixed = true,
                   description = "e_i: Measurement error on investment", tex_label = "e_i")

    # Setting steady-state parameters
    #nx = get_setting(m, :nx)
    nx = get_setting(m, :nx1_state) + get_setting(m, :nx2_state)
    ns = get_setting(m, :ns)

    # Steady state grids for functional/distributional variables and market-clearing discount rate
    m <= SteadyStateParameterGrid(:lstar, fill(NaN, nx*ns),
                                  description = "Steady-state expected discounted
                                  marginal utility of consumption", tex_label = "l_*")
    m <= SteadyStateParameterGrid(:cstar, fill(NaN, nx*ns),
                                  description = "Steady-state consumption", tex_label = "c_*")
    m <= SteadyStateParameterGrid(:μstar, fill(NaN, nx*ns), description = "Steady-state" *
                                  " cross-sectional density of cash on hand",
                                  tex_label = "\\mu_*")
    m <= SteadyStateParameter(:βstar, NaN, description = "Steady-state discount factor",
                              tex_label = "\\beta_*")
end

"""
```
aggregate_steadystate!(m::HetDSGEGovDebt)
```
Steady state of aggregate scalar variables (Break out the "analytic" steady-state solution here instead of in steadystate!() since some of these parameters are required to construct the grids that are used in the functional/distributional steady state variable calculation (where init_grids! is called before steadystate! in model construction).
"""
function aggregate_steadystate!(m::HetDSGEGovDebt)
    # FOR POSTERITY
    #=m <= SteadyStateParameter(:rprior, 100*(m[:r] - m[:γ]),
                              description= "100*(r* - γ)",
                              tex_label= "100*(r^{HetDSGE}-\\gamma^{FRBNY})")
    =#
    m <= SteadyStateParameter(:Rkstar, m[:r] + m[:δ],
                              description = "Rental rate on capital", tex_label = "Rk_*")
    m <= SteadyStateParameter(:ωstar,
                              (m[:α]^(m[:α]/(1-m[:α])))*(1-m[:α])*m[:Rkstar]^(-m[:α]/(1-m[:α])),
                              description = "Real wage", tex_label = "\\omega_*")
    m <= SteadyStateParameter(:klstar, (m[:α]/(1-m[:α]))*(m[:ωstar]/m[:Rkstar])*exp(m[:γ]),
                              description = "Capital/Labor ratio", tex_label = "kl_*")
    m <= SteadyStateParameter(:kstar, m[:klstar]*m[:H], description = "Capital",
                              tex_label = "k_*")
    m <= SteadyStateParameter(:xstar, (1-(1-m[:δ])*exp(-m[:γ]))*m[:kstar],
                              description = "Investment", tex_label = "x_*")
    m <= SteadyStateParameter(:ystar, (exp(-m[:α]*m[:γ])*(m[:kstar]^m[:α])*m[:H]^(1-m[:α])),
                              description = "GDP", tex_label = "y_*")
    m <= SteadyStateParameter(:bg, m[:BoverY]*m[:ystar]*exp(m[:γ]), description = "Govt Debt",
                              tex_label = "bg")
    m <= SteadyStateParameter(:Tg, (exp(-m[:γ]) - 1/(1+m[:r]))*m[:bg] + (1-(1/m[:g]))*m[:ystar],
                              description = "Net lump sump taxes", tex_label = "Tg")
    m <= SteadyStateParameter(:Tstar, m[:Rkstar]*m[:kstar]*exp(-m[:γ]) - m[:xstar] - m[:Tg],
                              description = "Net transfer to households", tex_label = "T_*")
    return m
end

"""
```
init_grids!(m::HetDSGEGovDebt)
```
"""
function init_grids!(m::HetDSGEGovDebt)


    nx      = get_setting(m, :nx)
    ns      = get_setting(m, :ns)
    λ       = get_setting(m, :λ)

    grids = OrderedDict()

    # Skill grid
    f, sgrid, swts, sscale = persistent_skill_process(m[:sH_over_sL].value, m[:pLH].value,
                                                      m[:pHL].value, get_setting(m, :ns))
    # Markov transition matrix for skill
    grids[:sgrid] = Grid(sgrid, swts, sscale)
    grids[:fgrid] = f

    xgrid, xwts, xlo, xhi, xscale = cash_grid(sgrid, m[:ωstar].value, m[:H].value,
                                              m[:r].scaledvalue, m[:η].value, m[:γ].scaledvalue,
                                              m[:Tstar].value, m[:zlo].value, nx)

    grids[:xgrid] = Grid(uniform_quadrature(xscale), xlo, xhi, nx, scale = xscale)

    m <= Setting(:xlo, xlo)
    m <= Setting(:xhi, xhi)
    m <= Setting(:xscale, xscale)

    # Total grid vectorized across both dimensions
    grids[:sgrid_total] = kron(sgrid, ones(nx))
    grids[:xgrid_total] = kron(ones(ns), grids[:xgrid].points)
    grids[:weights_total] = kron(swts, grids[:xgrid].weights)

    m.grids = grids
end

function model_settings!(m::HetDSGEGovDebt)
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

    m <= Setting(:policy_damp, 0.5, "Dampening parameter for policy function iteration")

    # Mollifier setting parameters
    m <= Setting(:In, 0.443993816237631, "Normalizing constant for the mollifier")

    # s: Skill Distribution/ "Units of effective labor" Grid Setup
    m <= Setting(:ns, 2, "Skill distribution grid points")
    m <= Setting(:λ, 2.0, "λ parameter in the Tauchen distribution calculation")

    # x: Cash on Hand Grid Setup
    m <= Setting(:nx1_state, 300, "Cash on hand distribution grid points (Lo)")
    m <= Setting(:nx2_state, 300, "Cash on hand distribution grid points (hi)")
    m <= Setting(:nx1_jump, 300, "Cash on hand distribution grid points (Lo)")
    m <= Setting(:nx2_jump, 300, "Cash on hand distribution grid points (hi)")
    m <= Setting(:nx, 300, "Cash on hand distribution grid points")

    m <= Setting(:binsize, 4) # Setting binsize=1 gives us what we had before doing the binning reduction
    m <= Setting(:poor_man_reduc, true) #note that we're actually doing more than the "poor man reduction" now however this turns ont both poorman truncation and binning reduction

    # Set targets
    m <= Setting(:targets, [0.16, 0.10],
                 "Targets for: [MPC, proportion of individuals with 0 income]")
    m <= Setting(:calibration_targets, [0.7, 0.23],
                 "Targets for: [var(log(annual income)), var(one year changes in " *
                 "log(annual income))]")
    m <= Setting(:calibration_targets_lb, [3.0, 1e-18], "Lower bounds on [sH/sL, zlo]")
    m <= Setting(:calibration_targets_ub, [9.0, 0.8-eps()], "Upper bounds on [sH/sL, zlo]")

    # Important settings for likelihood penalty
    m <= Setting(:use_likelihood_penalty, true)
    m <= Setting(:ψ_likelihood, 1.0,
                 "Multiplier on likelihood in penalty function")
    m <= Setting(:ψ_penalty, 1.0,
                 "Multiplier on likelihood in penalty function")

    m <= Setting(:target_vars, [:mpc, :pc0], "Symbols of variables we're targeting")
    m <= Setting(:target_σt, [0.2, 0.1], "Target \\sigma_t for MPC and pc0")

    m <= Setting(:steady_state_only, false, "Testing setting")
    m <= Setting(:auto_reject, false, "This flag is set when policy function doesn't converge")

    # Misc
    m <= Setting(:trunc_distr, false)
    m <= Setting(:rescale_weights, true)
    m <= Setting(:mindens, 1e-8)

    # Function-valued variables include distributional variables
    m <= Setting(:n_function_valued_backward_looking_states, 1, "Number of function-valued" *
                 " backward looking state variables")
    m <= Setting(:n_backward_looking_distributional_vars, 1,
                 "Number of state variables that are distributional variables.")
    m <= Setting(:n_function_valued_jumps, 1, "Number of function-valued jump variables")
    m <= Setting(:n_jump_distributional_vars, 1,
                 "Number of jump variables that are distributional variables.")

    # Note, these settings assume normalization.
    # The n degrees of freedom removed depends on the distributions/dimensions
    # of heterogeneity that we have discretized over, in this case,
    # cash on hand and the skill distribution. In general the rule of
    # thumb is, remove one degree of freedom for the first endogenous distribution (cash on
    # hand), then one additional degree of freedom for each exogenous distribution (skill
    # distribution). Multiple endogenous distributions only permit removing a single degree
    # of freedom since it is then non-trivial to obtain the marginal distributions.
    m <= Setting(:n_degrees_of_freedom_removed_state, 2,
                 "Number of degrees of freedom from the distributional variables to remove.")
    m <= Setting(:n_degrees_of_freedom_removed_jump, 0,
                 "Number of degrees of freedom from the distributional variables to remove.")
end

function setup_indices!(m::HetDSGEGovDebt)
    nx = get_setting(m, :nx)
    ns = get_setting(m, :ns)
    endo = m.endogenous_states_unnormalized
    eqconds = m.equilibrium_conditions
    nxns_state = (get_setting(m, :nx1_state) + get_setting(m, :nx2_state)) #*ns
    nxns_jump = (get_setting(m, :nx1_jump) + get_setting(m, :nx2_jump)) #*ns

    # Endogenous function-valued states
    endo[:kf′_t]    = 1:nxns_state    #  combination of lagged ell function and lagged m function that predicts m

    #endogenous scalar-valued states
    endo[:k′_t]   = nxns_state+1:nxns_state+1             # capital –dont get confused with steadystate object K
    endo[:R′_t1]  = nxns_state+2:nxns_state+2              # lagged real interest rate
    endo[:i′_t1]  = nxns_state+3:nxns_state+3             # lagged nominal interest rate
    endo[:y′_t1]  = nxns_state+4:nxns_state+4             # lagged gdp
    endo[:w′_t1]  = nxns_state+5:nxns_state+5           # lag real wages
    endo[:I′_t1]  = nxns_state+6:nxns_state+6           # lag investment–don't get this confused with i, the nominal interest rate
    endo[:bg′_t]  = nxns_state+7:nxns_state+7        # govt debt
    # exogenous scalar-valued states:
    endo[:b′_t]   = nxns_state+8:nxns_state+8        # discount factor shock
    endo[:g′_t]   = nxns_state+9:nxns_state+9        # govt spending
    endo[:z′_t]   = nxns_state+10:nxns_state+10        # tfp growth
    endo[:μ′_t]   = nxns_state+11:nxns_state+11        # investment shock
    endo[:λ_w′_t] = nxns_state+12:nxns_state+12        # wage markup
    endo[:λ_f′_t] = nxns_state+13:nxns_state+13        # price markup
    endo[:rm′_t]  = nxns_state+14:nxns_state+14        # monetary policy shock
    #endo[:c′_t1]  = nxns+14:nxns+14        # lagged consumption

    # function-valued jumps
    endo[:l′_t]  = nxns_state+15:nxns_state+nxns_jump+14 # ell function

    #scalar-valued jumps
    endo[:R′_t]   = nxns_state+nxns_jump+15:nxns_state+nxns_jump+15        # real interest rate
    endo[:i′_t]   = nxns_state+nxns_jump+16:nxns_state+nxns_jump+16        # nominal interest rate
    endo[:t′_t]   = nxns_state+nxns_jump+17:nxns_state+nxns_jump+17        # transfers + dividends
    endo[:w′_t]   = nxns_state+nxns_jump+18:nxns_state+nxns_jump+18        # real wage
    endo[:L′_t]   = nxns_state+nxns_jump+19:nxns_state+nxns_jump+19        # hours worked
    endo[:π′_t]   = nxns_state+nxns_jump+20:nxns_state+nxns_jump+20        # inflation
    endo[:π_w′_t] = nxns_state+nxns_jump+21:nxns_state+nxns_jump+21        # nominal wage inflation
    endo[:margutil′_t]  = nxns_state+nxns_jump+22:nxns_state+nxns_jump+22        # average marginal utility
    endo[:y′_t]   = nxns_state+nxns_jump+23:nxns_state+nxns_jump+23       # gdp
    endo[:I′_t]   = nxns_state+nxns_jump+24:nxns_state+nxns_jump+24        # investment
    endo[:mc′_t]  = nxns_state+nxns_jump+25:nxns_state+nxns_jump+25        # marginal cost - this is ζ in HetDSGEGovDebtₖd.pdf
    endo[:Q′_t]   = nxns_state+nxns_jump+26:nxns_state+nxns_jump+26        # Tobin's qfunction
    endo[:capreturn′_t] = nxns_state+nxns_jump+27:nxns_state+nxns_jump+27        # return on capital
    endo[:tg′_t] = nxns_state+nxns_jump+28:nxns_state+nxns_jump+28
    #endo[:c′_t] = 2*nxns+27:2*nxns+27        # consumption

    nvars = 2*nxns_jump+28
    nscalars = 28 # num eqs which output scalars
    nyscalars = 14 # num scalar jumps
    nxscalars = nscalars - nyscalars # num scalar states

    # create objects needed for solve.jl
    # we will order function blocks as follows:
    # 1. all function blocks which output a function (first real eqs, then lags)
    # 2. all function blocks which map functions to scalars
    # 3. all scalar blocks involving endogenous vbls (first real eqs, then lags)
    # 4. shock processes
    funops = 1:2 # which operators output a function

    # function blocks which output a function
    eqconds[:eq_euler]              = 1:nxns_state
    eqconds[:eq_kolmogorov_fwd]     = nxns_state+1:2*nxns_state

    # function blocks which map functions to scalars
    eqconds[:eq_market_clearing]    = 2*nxns_state+1:2*nxns_state+1
    eqconds[:eq_lambda]             = 2*nxns_state+2:2*nxns_state+2
    #scalar blocks involving endogenous variables
    eqconds[:eq_transfers]            = 2*nxns_state+3:2*nxns_state+3 # transfers
    eqconds[:eq_investment]           = 2*nxns_state+4:2*nxns_state+4 # investment
    eqconds[:eq_tobin_q]              = 2*nxns_state+5:2*nxns_state+5 # tobin's q
    eqconds[:eq_capital_accumulation] = 2*nxns_state+6:2*nxns_state+6 # capital accumulation
    eqconds[:eq_wage_phillips] = 2*nxns_state+7:2*nxns_state+7 # wage phillips curve
    eqconds[:eq_price_phillips] = 2*nxns_state+8:2*nxns_state+8 # price phillips curve
    eqconds[:eq_marginal_cost]  = 2*nxns_state+9:2*nxns_state+9 # marginal cost
    eqconds[:eq_gdp]  = 2*nxns_state+10:2*nxns_state+10 # gdp
    eqconds[:eq_optimal_kl] = 2*nxns_state+11:2*nxns_state+11 # optimal K/L ratio
    eqconds[:eq_taylor] = 2*nxns_state+12:2*nxns_state+12 # taylor rule
    eqconds[:eq_fisher] = 2*nxns_state+13:2*nxns_state+13 # fisher eqn
    eqconds[:eq_nominal_wage_inflation] = 2*nxns_state+14:2*nxns_state+14 # nominal wage inflation
    eqconds[:eq_fiscal_rule] = 2*nxns_state+15:2*nxns_state+15
    eqconds[:eq_g_budget_constraint] = 2*nxns_state+16:2*nxns_state+16

    # lagged variables
    eqconds[:LR] = 2*nxns_state+17:2*nxns_state+17 # LR
    eqconds[:LI] = 2*nxns_state+18:2*nxns_state+18 # LI
    eqconds[:LY] = 2*nxns_state+19:2*nxns_state+19 # LY
    eqconds[:LW]  = 2*nxns_state+20:2*nxns_state+20 # LW
    eqconds[:LX] = 2*nxns_state+21:2*nxns_state+21 # LX
    # shocks
    eqconds[:eq_b] = 2*nxns_state+22:2*nxns_state+22 # discount factor B
    eqconds[:eq_g] = 2*nxns_state+23:2*nxns_state+23 # govt spending G
    eqconds[:eq_z] = 2*nxns_state+24:2*nxns_state+24 # tfp growth Z
    eqconds[:eq_μ] = 2*nxns_state+25:2*nxns_state+25 # investment MU
    eqconds[:eq_λ_w] = 2*nxns_state+26:2*nxns_state+26 # wage mkup LAMW
    eqconds[:eq_λ_f] = 2*nxns_state+27:2*nxns_state+27 # price mkup LAMF
    eqconds[:eq_rm] = 2*nxns_state+28:2*nxns_state+28 # monetary policy MON
    #eqconds[:eq_consumption] = 2*nxns+27:2*nxns+27 # monetary policy MON

    # Total grid x*s
    m <= Setting(:n_state, (get_setting(m, :nx1_state) +get_setting(m, :nx2_state)),
                 "Total grid size, multiplying across grid dimensions.")
m <= Setting(:n_jump, (get_setting(m, :nx1_jump) +get_setting(m, :nx2_jump)),
             "Total grid size, multiplying across grid dimensions.")
    m <= Setting(:nvars, 2*get_setting(m, :n_state) + 28, "num variables")
    m <= Setting(:nscalars, 28, " # num eqs which output scalars")
    m <= Setting(:nyscalars, 14, "num scalar jumps")
    m <= Setting(:nxscalars, get_setting(m, :nscalars) - get_setting(m, :nyscalars),
                 "num scalar states")
    m.endogenous_states = deepcopy(endo)
end

function init_states_and_jumps!(m::AbstractModel, states::Vector{Symbol}, jumps::Vector{Symbol}, states_only::Bool = false)
    endo = m.endogenous_states_unnormalized

    m <= Setting(:states, states)
    m <= Setting(:jumps, jumps)

    state_indices = stack_indices(endo, states)
    jump_indices = stack_indices(endo, jumps)

    m <= Setting(:state_indices, state_indices, "Indices of m.endogenous_states " *
                 "corresponding to backward-looking state variables")
    m <= Setting(:jump_indices, jump_indices, "Indices of m.endogenous_states " *
                 "corresponding to jump
                 variables")

    n_dof_removed_state = get_setting(m, :n_degrees_of_freedom_removed_state)
    n_dof_removed_jump  = get_setting(m, :n_degrees_of_freedom_removed_jump)


    ####################################################
    # Calculating the number of backward-looking states
    ####################################################
    n_backward_looking_distr_vars = get_setting(m, :n_backward_looking_distributional_vars)
    m <= Setting(:backward_looking_states_normalization_factor,
                 n_dof_removed_state*n_backward_looking_distr_vars,
                 "Number of dimensions removed from the backward looking state variables
                  for the normalization.")

    nxns_state = (get_setting(m, :nx1_state) + get_setting(m, :nx2_state))*get_setting(m, :ns)#get_setting(m, :n) #nx)*get_setting(m, :ns)
    nxns_jump = (get_setting(m, :nx1_jump) + get_setting(m, :nx2_jump))*get_setting(m, :ns)#get_setting(m, :n) #nx)*get_setting(m, :ns)
    nx1_state = get_setting(m, :nx1_state)
    nx2_state = get_setting(m, :nx2_state)
    nx1_jump = get_setting(m, :nx1_jump)
    nx2_jump = get_setting(m, :nx2_jump)
    n_backward_looking_vars = length(get_setting(m, :state_indices))
    n_backward_looking_function_valued_vars = get_setting(m,
                                               :n_function_valued_backward_looking_states)
    n_backward_looking_scalar_vars = get_setting(m, :nxscalars) #=n_backward_looking_vars -
        (nx1 + nx2)*n_backward_looking_function_valued_vars =# #n_backward_looking_function_valued_vars

    m <= Setting(:n_backward_looking_states, (nx1_state + nx2_state)*n_backward_looking_distr_vars +
                 n_backward_looking_scalar_vars -
                 get_setting(m, :backward_looking_states_normalization_factor),
                 "Number of state variables, in the true sense (fully
                  backward looking) accounting for the discretization across the grid
                  of function-valued variables and the normalization of
                  distributional variables.")

    ##################################
    # Calculating the number of jumps
    ##################################
    n_jump_distr_vars = get_setting(m, :n_jump_distributional_vars)
    m <= Setting(:jumps_normalization_factor,
                 n_dof_removed_jump*n_jump_distr_vars, "number of dimensions removed from" *
                 " the jump variables for the normalization.")

    n_jump_vars = length(get_setting(m, :jump_indices))
    n_jump_function_valued_vars = get_setting(m, :n_function_valued_jumps)
    n_jump_scalar_vars = get_setting(m, :nyscalars) #n_jump_vars - (nx1+nx2)*n_jump_function_valued_vars

    m <= Setting(:n_jumps, (nx1_jump+nx2_jump)*n_jump_function_valued_vars +
                 n_jump_scalar_vars - get_setting(m, :jumps_normalization_factor),
                 "Number of jump variables (forward looking) accounting for
                  the discretization across the grid of function-valued variables and the
                  normalization of distributional variables.")

    m <= Setting(:n_model_states, get_setting(m, :n_backward_looking_states) +
                 get_setting(m, :n_jumps),
                 "Number of 'states' in the state space model. Because backward and forward
                 looking variables need to be explicitly tracked for the Klein solution
                 method, we have n_states and n_jumps")

end

function reset_grids!(m)
    m <= Setting(:nx1_state, 300)
    m <= Setting(:nx2_state, 300)
    m <= Setting(:nx1_jump, 300)
    m <= Setting(:nx2_jump, 300)


    setup_indices!(m)

    init_states_and_jumps!(m, get_setting(m, :states), get_setting(m, :jumps))
    init_grids!(m)

    # So that the indices of m.endogenous_states reflect the normalization
    normalize_model_state_indices!(m)

    endogenous_states_augmented = [:i_t1, :c_t, :c_t1]
    for (i,k) in enumerate(endogenous_states_augmented)
        m.endogenous_states_augmented[k] = i + first(collect(values(m.endogenous_states))[end])
    end
    m <= Setting(:n_model_states_augmented, get_setting(m, :n_model_states) +
                 length(m.endogenous_states_augmented))

end
