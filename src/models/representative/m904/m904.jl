"""
```
Model904{T} <: AbstractDSGEModel{T}
```

The `Model904` type defines the structure of the New York Fed DSGE model.

### Fields

#### Parameters and Steady-States
* `parameters::Vector{AbstractParameter}`: Vector of all time-invariant model parameters.

* `steady_state::Vector{AbstractParameter}`: Model steady-state values, computed as a function of elements of
  `parameters`.

* `keys::OrderedDict{Symbol,Int}`: Maps human-readable names for all model parameters and
  steady-states to their indices in `parameters` and `steady_state`.

#### Inputs to Measurement and Equilibrium Condition Equations

The following fields are dictionaries that map human-readible names to row and column
indices in the matrix representations of of the measurement equation and equilibrium
conditions.

* `endogenous_states::OrderedDict{Symbol,Int}`: Maps each state to a column in the measurement and
  equilibrium condition matrices.

* `exogenous_shocks::OrderedDict{Symbol,Int}`: Maps each shock to a column in the measurement and
  equilibrium condition matrices.

* `expected_shocks::OrderedDict{Symbol,Int}`: Maps each expected shock to a column in the
  measurement and equilibrium condition matrices.

* `equilibrium_conditions::OrderedDict{Symbol,Int}`: Maps each equlibrium condition to a row in the
  model's equilibrium condition matrices.

* `endogenous_states_augmented::OrderedDict{Symbol,Int}`: Maps lagged states to their columns in
  the measurement and equilibrium condition equations. These are added after Gensys solves the
  model.

* `observables::OrderedDict{Symbol,Int}`: Maps each observable to a row in the model's measurement
  equation matrices.

* `pseudo_observables::OrderedDict{Symbol,Int}`: Maps each pseudo-observable to
  a row in the model's pseudo-measurement equation matrices.

#### Model Specifications and Settings

* `spec::String`: The model specification identifier, \"m904\", cached here for
  filepath computation.

* `subspec::String`: The model subspecification number, indicating that some
  parameters from the original model spec (\"ss0\") are initialized differently. Cached here for
  filepath computation.

* `settings::Dict{Symbol,Setting}`: Settings/flags that affect computation without changing
  the economic or mathematical setup of the model.

* `test_settings::Dict{Symbol,Setting}`: Settings/flags for testing mode

#### Other Fields

* `rng::MersenneTwister`: Random number generator. Can be is seeded to ensure
  reproducibility in algorithms that involve randomness (such as Metropolis-Hastings).

* `testing::Bool`: Indicates whether the model is in testing mode. If `true`, settings from
  `m.test_settings` are used in place of those in `m.settings`.

* `observable_mappings::OrderedDict{Symbol,Observable}`: A dictionary that
  stores data sources, series mnemonics, and transformations to/from model units.
  DSGE.jl will fetch data from the Federal Reserve Bank of
  St. Louis's FRED database; all other data must be downloaded by the
  user. See `load_data` and `Observable` for further details.

* `pseudo_observable_mappings::OrderedDict{Symbol,PseudoObservable}`: A
  dictionary that stores names and transformations to/from model units. See
  `PseudoObservable` for further details.
"""
mutable struct Model904{T} <: AbstractRepModel{T}
    parameters::ParameterVector{T}                       # vector of all time-invariant model parameters
    steady_state::ParameterVector{T}                     # model steady-state values
    keys::OrderedDict{Symbol,Int}                        # human-readable names for all the model
                                                         # parameters and steady-n_states

    endogenous_states::OrderedDict{Symbol,Int}           # these fields used to create matrices in the
    exogenous_shocks::OrderedDict{Symbol,Int}            # measurement and equilibrium condition equations.
    expected_shocks::OrderedDict{Symbol,Int}             #
    equilibrium_conditions::OrderedDict{Symbol,Int}      #
    endogenous_states_augmented::OrderedDict{Symbol,Int} #
    observables::OrderedDict{Symbol,Int}                 #
    pseudo_observables::OrderedDict{Symbol,Int}          #
    spec::String                                         # Model specification number (eg "m1006")
    subspec::String                                      # Model subspecification (eg "ss0")
    settings::Dict{Symbol,Setting}                       # Settings/flags for computation
    test_settings::Dict{Symbol,Setting}                  # Settings/flags for testing mode
    rng::MersenneTwister                                 # Random number generator
    testing::Bool                                        # Whether we are in testing mode or not

    observable_mappings::OrderedDict{Symbol, Observable}
    pseudo_observable_mappings::OrderedDict{Symbol, PseudoObservable}
end

description(m::Model904) = "New York Fed DSGE Model m904, $(m.subspec)"

"""
`init_model_indices!(m::Model904)`

Arguments:
`m:: Model904`: a model object

Description:
Initializes indices for all of `m`'s states, shocks, and equilibrium conditions.
"""
function init_model_indices!(m::Model904)
    # Endogenous states
    endogenous_states = [[
        :y_t, :c_t, :i_t, :qk_t, :k_t, :kbar_t, :u_t, :rk_t, :Rtil_k_t, :n_t, :mc_t,
        :π_t, :μ_ω_t, :w_t, :L_t, :R_t, :g_t, :b_t, :μ_t, :z_t, :λ_f_t, :λ_f_t1,
        :λ_w_t, :λ_w_t1, :rm_t, :σ_ω_t, :μ_e_t, :γ_t, :π_star_t, :Ec_t, :Eqk_t, :Ei_t,
        :Eπ_t, :EL_t, :Erk_t, :Ew_t, :ERtil_k_t, :y_f_t, :c_f_t, :i_f_t, :qk_f_t, :k_f_t,
        :kbar_f_t, :u_f_t, :rk_f_t, :w_f_t, :L_f_t, :r_f_t, :Ec_f_t, :Eqk_f_t, :Ei_f_t,
        :EL_f_t, :Erk_f_t, :ztil_t];
        [Symbol("rm_tl$i") for i = 1:n_anticipated_shocks(m)]]

    # Exogenous shocks
    exogenous_shocks = [[
        :g_sh, :b_sh, :μ_sh, :z_sh, :λ_f_sh, :λ_w_sh, :rm_sh, :σ_ω_sh, :μ_e_sh,
        :γ_sh, :π_star_sh];
        [Symbol("rm_shl$i") for i = 1:n_anticipated_shocks(m)]]

    # Expectations shocks
    expected_shocks = [
        :Ec_sh, :Eqk_sh, :Ei_sh, :Eπ_sh, :EL_sh, :Erk_sh, :Ew_sh, :ERktil_sh, :Ec_f_sh,
        :Eqk_f_sh, :Ei_f_sh, :EL_f_sh, :Erk_f_sh]

    # Equilibrium conditions
    equilibrium_conditions = [[
        :eq_euler, :eq_inv, :eq_capval, :eq_spread, :eq_nevol, :eq_output, :eq_caputl, :eq_capsrv, :eq_capev,
        :eq_mkupp, :eq_phlps, :eq_caprnt, :eq_msub, :eq_wage, :eq_mp, :eq_res, :eq_g, :eq_b, :eq_μ, :eq_z,
        :eq_λ_f, :eq_λ_w, :eq_rm, :eq_σ_ω, :eq_μ_e, :eq_γ, :eq_λ_f1, :eq_λ_w1, :eq_Ec,
        :eq_Eqk, :eq_Ei, :eq_Eπ, :eq_EL, :eq_Erk, :eq_Ew, :eq_ERktil, :eq_euler_f, :eq_inv_f,
        :eq_capval_f, :eq_output_f, :eq_caputl_f, :eq_capsrv_f, :eq_capev_f, :eq_mkupp_f,
        :eq_caprnt_f, :eq_msub_f, :eq_res_f, :eq_Ec_f, :eq_Eqk_f, :eq_Ei_f, :eq_EL_f, :eq_Erk_f,
        :eq_ztil, :eq_π_star];
        [Symbol("eq_rml$i") for i=1:n_anticipated_shocks(m)]]

    # Additional states added after solving model
    # Lagged states and observables measurement error
    endogenous_states_augmented = [
        :y_t1, :c_t1, :i_t1, :w_t1, :π_t1_dup, :L_t1, :Et_π_t]

    # Observables
    observables = keys(m.observable_mappings)

    # Pseudo-observables
    pseudo_observables = keys(m.pseudo_observable_mappings)

    for (i,k) in enumerate(endogenous_states);           m.endogenous_states[k]           = i end
    for (i,k) in enumerate(exogenous_shocks);            m.exogenous_shocks[k]            = i end
    for (i,k) in enumerate(expected_shocks);             m.expected_shocks[k]             = i end
    for (i,k) in enumerate(equilibrium_conditions);      m.equilibrium_conditions[k]      = i end
    for (i,k) in enumerate(endogenous_states);           m.endogenous_states[k]           = i end
    for (i,k) in enumerate(endogenous_states_augmented); m.endogenous_states_augmented[k] = i+length(endogenous_states) end
    for (i,k) in enumerate(observables);                 m.observables[k]                 = i end
    for (i,k) in enumerate(pseudo_observables);          m.pseudo_observables[k]          = i end
end


function Model904(subspec::String="ss9";
                  custom_settings::Dict{Symbol, Setting} = Dict{Symbol, Setting}(),
                  testing = false)

    # Model-specific specifications
    spec               = split(basename(@__FILE__),'.')[1]
    subspec            = subspec
    settings           = Dict{Symbol,Setting}()
    test_settings      = Dict{Symbol,Setting}()
    rng                = MersenneTwister(0)

    # initialize empty model
    m = Model904{Float64}(
            # model parameters and steady state values
            Vector{AbstractParameter{Float64}}(), Vector{Float64}(), OrderedDict{Symbol,Int}(),

            # model indices
            OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(),

            spec,
            subspec,
            settings,
            test_settings,
            rng,
            testing,
            OrderedDict{Symbol,Observable}(),
            OrderedDict{Symbol,PseudoObservable}())

    # Set settings
    model_settings!(m)
    default_test_settings!(m)
    for custom_setting in values(custom_settings)
        m <= custom_setting
    end

    # Set observable and pseudo-observable transformations
    init_observable_mappings!(m)
    init_pseudo_observable_mappings!(m)

    # Initialize parameters
    init_parameters!(m)

    init_model_indices!(m)
    init_subspec!(m)
    steadystate!(m)

    return m
end

"""
```
init_parameters!(m::Model904)
```

Initializes the model's parameters, as well as empty values for the steady-state
parameters (in preparation for `steadystate!(m)` being called to initialize
those).
"""
function init_parameters!(m::Model904)
    m <= parameter(:α,      0.1687, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     Normal(0.30, 0.05),         fixed=false,
                   description="α: Capital elasticity in the intermediate goods sector's production function (also known as the capital share).",
                   tex_label="\\alpha")

    m <= parameter(:ζ_p,   0.7467, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     BetaAlt(0.5, 0.1),          fixed=false,
                   description="ζ_p: The Calvo parameter. In every period, intermediate goods producers optimize prices with probability (1-ζ_p). With probability ζ_p, prices are adjusted according to a weighted average of the previous period's inflation (π_t1) and steady-state inflation (π_star).",
                   tex_label="\\zeta_p")

    m <= parameter(:ι_p,   0.2684, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     BetaAlt(0.5, 0.15),         fixed=false,
                   description="ι_p: The weight attributed to last period's inflation in price indexation. (1-ι_p) is the weight attributed to steady-state inflation.",
                   tex_label="\\iota_p")

    m <= parameter(:δ,      0.025,  fixed=true,
                   description="δ: The capital depreciation rate.", tex_label="\\delta" )

    m <= parameter(:Upsilon,  1.000,  (0., 10.),     (1e-5, 0.),      ModelConstructors.Exponential(),    GammaAlt(1., 0.5),          fixed=true,
                   description="Υ: The trend evolution of the price of investment goods relative to consumption goods. Set equal to 1.",
                   tex_label="\\mathcal{\\Upsilon}")

    m <= parameter(:Φ,   1.4933, (1., 10.),     (1.00, 10.00),   ModelConstructors.Exponential(),    Normal(1.25, 0.12),         fixed=false,
                   description="Φ: Fixed costs.",
                   tex_label="\\Phi")

    m <= parameter(:S′′,       3.3543, (-15., 15.),   (-15., 15.),     Untransformed(),  Normal(4., 1.5),            fixed=false,
                   description="S'': The second derivative of households' cost of adjusting investment.", tex_label="S\\prime\\prime")

    m <= parameter(:h,        0.4656, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     BetaAlt(0.7, 0.1),          fixed=false,
                   description="h: Consumption habit persistence.", tex_label="h")

    m <= parameter(:ppsi,     0.7614, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     BetaAlt(0.5, 0.15),         fixed=false,
                   description="ppsi: Utilization costs.", tex_label="\\psi")

    m <= parameter(:ν_l,     1.0647, (1e-5, 10.),   (1e-5, 10.),     ModelConstructors.Exponential(),    Normal(2, 0.75),            fixed=false,
                   description="ν_l: The coefficient of relative risk aversion on the labor
                   term of households' utility function.", tex_label="\\nu_l")

    m <= parameter(:ζ_w,   0.7922, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     BetaAlt(0.5, 0.1),          fixed=false,
                   description="ζ_w: (1-ζ_w) is the probability with which households can freely choose wages in each period. With probability ζ_w, wages increase at a geometrically weighted average of the steady state rate of wage increases and last period's productivity times last period's inflation.",
                   tex_label="\\zeta_w")

    m <= parameter(:ι_w,   0.5729, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     BetaAlt(0.5, 0.15),         fixed=false,
                   description="ι_w: No description available.",
                   tex_label="\\iota_w")

    m <= parameter(:λ_w,      1.5000,                                                                               fixed=true,
                   description="λ_w: The wage markup, which affects the elasticity of substitution between differentiated labor services.",
                   tex_label="\\lambda_w")

    m <= parameter(:β,      0.7420, (1e-5, 10.),   (1e-5, 10.),     ModelConstructors.Exponential(),    GammaAlt(0.25, 0.1),        fixed=false,  scaling = x -> 1/(1 + x/100),
                   description="β: Discount rate.",
                   tex_label="100(\\beta^{-1} - 1)")

    m <= parameter(:ψ1,     1.8678, (1e-5, 10.),   (1e-5, 10.00),   ModelConstructors.Exponential(),    Normal(1.5, 0.25),          fixed=false,
                   description="ψ₁: Weight on inflation gap in monetary policy rule.",
                   tex_label="\\psi_1")

    m <= parameter(:ψ2,     0.0715, (-0.5, 0.5),   (-0.5, 0.5),     Untransformed(),  Normal(0.12, 0.05),         fixed=false,
                   description="ψ₂: Weight on output gap in monetary policy rule.",
                   tex_label="\\psi_2")

    m <= parameter(:ψ3,     0.2131, (-0.5, 0.5),   (-0.5, 0.5),     Untransformed(),  Normal(0.12, 0.05),         fixed=false,
                   description="ψ₃: Weight on rate of change of output gap in the monetary policy rule.",
                   tex_label="\\psi_3")

    m <= parameter(:π_star,   0.6231, (1e-5, 10.),   (1e-5, 10.),     ModelConstructors.Exponential(),    GammaAlt(0.75, 0.4),        fixed=false,  scaling = x -> 1 + x/100,
                   description="π_star: The steady-state rate of inflation.",
                   tex_label="\\pi_*")

    m <= parameter(:σ_c,   1.5073, (1e-5, 10.),   (1e-5, 10.),     ModelConstructors.Exponential(),    Normal(1.5, 0.37),          fixed=false,
                   description="σ_c: No description available.",
                   tex_label="\\sigma_{c}")

    m <= parameter(:ρ,      0.8519, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     BetaAlt(0.75, 0.10),        fixed=false,
                   description="ρ: The degree of inertia in the monetary policy rule.",
                   tex_label="\\rho")

    m <= parameter(:ϵ_p,     10.000,                                                                               fixed=true,
                   description="ϵ_p: No description available.",
                   tex_label="\\varepsilon_{p}")

    m <= parameter(:ϵ_w,     10.000,                                                                               fixed=true,
                   description="ϵ_w: No description available.",
                   tex_label="\\varepsilon_{w}")

    # financial frictions parameters
    m <= parameter(:Fω,      0.0300, (1e-5, 0.99999), (1e-5, 0.99),  SquareRoot(),    BetaAlt(0.03, 0.01),         fixed=true,    scaling = x -> 1 - (1-x)^0.25,
                   description="F(ω): The cumulative distribution function of ω (idiosyncratic iid shock that increases or decreases entrepreneurs' capital).",
                   tex_label="F(\\omega)")

    m <= parameter(:spr,     2.0, (0., 100.),      (1e-5, 0.),    ModelConstructors.Exponential(),   GammaAlt(2., 0.1),           fixed=false,  scaling = x -> (1 + x/100)^0.25,
                   description="spr_*: Steady-state spreads.",
                   tex_label="spr_*")

    m <= parameter(:ζ_spb, 0.05, (1e-5, 0.99999), (1e-5, 0.99),  SquareRoot(),    BetaAlt(0.05, 0.005),        fixed=false,
                   description="ζ_spb: The elasticity of the expected exess return on capital (or 'spread') with respect to leverage.",
                   tex_label="\\zeta_{spb}")

    m <= parameter(:γ_star, 0.9900, (1e-5, 0.99999), (1e-5, 0.99),  SquareRoot(),    BetaAlt(0.99, 0.002),        fixed=true,
                   description="γ_star: No description available.",
                   tex_label="\\gamma_*")

    # exogenous processes - level
    m <= parameter(:γ,      0.3085, (-5.0, 5.0),     (-5., 5.),     Untransformed(), Normal(0.4, 0.1),            fixed=false, scaling = x -> x/100,
                   description="γ: The log of the steady-state growth rate of technology.",
                   tex_label="100\\gamma")

    m <= parameter(:Lmean,  -1.0411, (-1000., 1000.), (-1e3, 1e3),   Untransformed(), Normal(-45, 5),              fixed=false,
                   description="Lmean: No description available.",
                   tex_label="\\bar{L}")

    m <= parameter(:g_star,    0.1800,                                                                               fixed=true,
                   description="g_star: No description available.",
                   tex_label="g_*")

    # exogenous processes - autocorrelation
    m <= parameter(:ρ_g,      0.9930, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(),    BetaAlt(0.5, 0.2),           fixed=false,
                   description="ρ_g: AR(1) coefficient in the government spending process.",
                   tex_label="\\rho_g")

    m <= parameter(:ρ_b,      0.8100, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(),    BetaAlt(0.5, 0.2),           fixed=false,
                   description="ρ_b: AR(1) coefficient in the intertemporal preference shifter process.",
                   tex_label="\\rho_b")

    m <= parameter(:ρ_μ,     0.7790, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(),    BetaAlt(0.5, 0.2),           fixed=false,
                   description="ρ_μ: AR(1) coefficient in capital adjustment cost process.",
                   tex_label="\\rho_{\\mu}")

    m <= parameter(:ρ_z,      0.9676, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(),    BetaAlt(0.5, 0.2),           fixed=false,
                   description="ρ_z: AR(1) coefficient in the technology process.",
                   tex_label="\\rho_z")

    m <= parameter(:ρ_λ_f,    0.8372, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(),    BetaAlt(0.5, 0.2),           fixed=false,
                   description="ρ_λ_f: AR(1) coefficient in the price mark-up shock process.",
                   tex_label="\\rho_{\\lambda_f}")

    m <= parameter(:ρ_λ_w,    0.9853, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(),    BetaAlt(0.5, 0.2),           fixed=false,
                   description="ρ_λ_w: AR(1) coefficient in the wage mark-up shock process.",
                   tex_label="\\rho_{\\lambda_w}")

    #monetary policy shock - see eqcond
    m <= parameter(:ρ_rm,     0.9000, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(),    BetaAlt(0.5, 0.2),           fixed=false,
                   description="ρ_rm: AR(1) coefficient in the monetary policy shock process.",
                   tex_label="\\rho_{rm}")

    m <= parameter(:ρ_σ_w,   0.9, (1e-5, 0.9999), (1e-5, 0.9999),  SquareRoot(),    BetaAlt(0.75, 0.15),         fixed=false,
                   description="ρ_σ_w: The standard deviation of entrepreneurs' capital productivity follows an exogenous process with mean ρ_σ_w. Innovations to the process are called _spread shocks_.",
                   tex_label="\\rho_{\\sigma_\\omega}")

    m <= parameter(:ρ_μ_e,    0.7500, (1e-5, 0.99999), (1e-5, 0.99999),  SquareRoot(),    BetaAlt(0.75, 0.15),         fixed=true,
                   description="ρ_μ_e: Verification costs are a fraction μ_e of the amount the bank extracts from an entrepreneur in case of bankruptcy???? This doesn't seem right because μ_e isn't a process (p12 of PDF)",
                   tex_label="\\rho_{\\mu_e}")

    m <= parameter(:ρ_γ,   0.7500, (1e-5, 0.99999), (1e-5, 0.99),  SquareRoot(), BetaAlt(0.75, 0.15),         fixed=true,  description="ρ_γ: Autocorrelation coefficient on the γ shock.",              tex_label="\\rho_{\\gamma}")
    m <= parameter(:ρ_π_star,   0.9900, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(), BetaAlt(0.5, 0.2),           fixed=true,  description="ρ_π_star: Autocorrelation coefficient on the π_star shock.",  tex_label="\\rho_{\\pi^*}")

    # exogenous processes - standard deviation
    m <= parameter(:σ_g,      0.6090, (1e-8, 5.),      (1e-8, 5.),    ModelConstructors.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
                   description="σ_g: The standard deviation of the government spending process.",
                   tex_label="\\sigma_{g}")

    m <= parameter(:σ_b,      0.1818, (1e-8, 5.),      (1e-8, 5.),    ModelConstructors.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
                   description="σ_b: No description available.",
                   tex_label="\\sigma_{b}")

    m <= parameter(:σ_μ,     0.4601, (1e-8, 5.),      (1e-8, 5.),    ModelConstructors.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
                   description="σ_μ: The standard deviation of the exogenous marginal efficiency of investment shock process.",
                   tex_label="\\sigma_{\\mu}")

    m <= parameter(:σ_z,      0.4618, (1e-8, 5.),      (1e-8, 5.),    ModelConstructors.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
                   description="σ_z: No description available.",
                   tex_label="\\sigma_{z}")

    m <= parameter(:σ_λ_f,    0.1455, (1e-8, 5.),      (1e-8, 5.),    ModelConstructors.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
                   description="σ_λ_f: The mean of the process that generates the price elasticity of the composite good.  Specifically, the elasticity is (1+λ_{f,t})/(λ_{f_t}).",
                   tex_label="\\sigma_{\\lambda_f}")

    m <= parameter(:σ_λ_w,    0.2089, (1e-8, 5.),      (1e-8, 5.),    ModelConstructors.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
                   description="σ_λ_w: No description available.",
                   tex_label="\\sigma_{\\lambda_w}")

    m <= parameter(:σ_r_m,     0.2397, (1e-8, 5.),      (1e-8, 5.),    ModelConstructors.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
                   description="σ_r_m: No description available.",
                   tex_label="\\sigma_{rm}")

    m <= parameter(:σ_σ_ω,   0.2/4, (1e-7,100.),     (1e-5, 0.),    ModelConstructors.Exponential(),   RootInverseGamma(4., 0.05),  fixed=false,
                   description="σ_σ_ω: The standard deviation of entrepreneurs' capital productivity follows an exogenous process with standard deviation σ_σ_ω.",
                   tex_label="\\sigma_{\\sigma_\\omega}")

    m <= parameter(:σ_μ_e,    0.0000, (1e-7,100.),     (1e-5, 0.),    ModelConstructors.Exponential(), RootInverseGamma(4., 0.05),  fixed=true,  description="σ_μ_e: Standard deviation of the μ_e shock.",         tex_label="\\sigma_{\\mu_e}")
    m <= parameter(:σ_γ,   0.0000, (1e-7,100.),     (1e-5, 0.),    ModelConstructors.Exponential(),   RootInverseGamma(4., 0.01),  fixed=true,  description="σ_γ: Standard deviation of the γ shock.",              tex_label="\\sigma_{\\gamma}")
    m <= parameter(:σ_π_star,   0.03, (1e-8, 5.),      (1e-8, 5.),    ModelConstructors.Exponential(), RootInverseGamma(6., 0.03),  fixed=false, description="σ_π_star: Standard deviation of the π_star shock.", tex_label="\\sigma_{\\pi^*}")

    m <= parameter(:η_gz,       0.5632, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(),
                   BetaAlt(0.50, 0.20), fixed=false,
                   description="η_gz: Correlate g and z shocks.",
                   tex_label="\\eta_{gz}")

    m <= parameter(:η_λ_f,      0.7652, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(),
                   BetaAlt(0.50, 0.20),         fixed=false,
                   description="η_λ_f: Moving average component in the price markup shock.",
                   tex_label="\\eta_{\\lambda_f}")

    m <= parameter(:η_λ_w,      0.8936, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(),
                   BetaAlt(0.50, 0.20),         fixed=false,
                   description="η_λ_w: Moving average component in the wage markup shock.",
                   tex_label="\\eta_{\\lambda_w}")

    # steady states
    m <= SteadyStateParameter(:z_star,  NaN, description="No description available.", tex_label="\\z_*")
    m <= SteadyStateParameter(:rstar,   NaN, description="No description available.", tex_label="\\r_*")
    m <= SteadyStateParameter(:Rstarn,  NaN, description="No description available.", tex_label="\\R_*_n")
    m <= SteadyStateParameter(:r_k_star,  NaN, description="No description available.", tex_label="\\r^k_*")
    m <= SteadyStateParameter(:wstar,   NaN, description="No description available.", tex_label="\\w_*")
    m <= SteadyStateParameter(:Lstar,   NaN, description="No description available.", tex_label="\\L_*")
    m <= SteadyStateParameter(:kstar,   NaN, description="Effective capital that households rent to firms in the steady state.", tex_label="\\k_*")
    m <= SteadyStateParameter(:kbarstar, NaN, description="Total capital owned by households in the steady state.", tex_label="\\bar{k}_*")
    m <= SteadyStateParameter(:istar,  NaN, description="Detrended steady-state investment", tex_label="\\i_*")
    m <= SteadyStateParameter(:ystar,  NaN, description="No description available.", tex_label="\\y_*")
    m <= SteadyStateParameter(:cstar,  NaN, description="No description available.", tex_label="\\c_*")
    m <= SteadyStateParameter(:wl_c,   NaN, description="No description available.", tex_label="\\wl_c")
    m <= SteadyStateParameter(:nstar,  NaN, description="No description available.", tex_label="\\n_*")
    m <= SteadyStateParameter(:vstar,  NaN, description="No description available.", tex_label="\\v_*")
    m <= SteadyStateParameter(:ζ_spσ_ω,  NaN, description="No description available.", tex_label="\\zeta_{sp_{\\sigma_\\omega}}")
    m <= SteadyStateParameter(:ζ_spμ_e,   NaN, description="No description available.", tex_label="\\zeta_{sp_{\\mu_e}}")
    m <= SteadyStateParameter(:ζ_nRk,     NaN, description="No description available.", tex_label="\\zeta_{n_R_k}")
    m <= SteadyStateParameter(:ζ_nR,      NaN, description="No description available.", tex_label="\\zeta_{n_R}")
    m <= SteadyStateParameter(:ζ_nqk,     NaN, description="No description available.", tex_label="\\zeta_{n_q_k}")
    m <= SteadyStateParameter(:ζ_nn,      NaN, description="No description available.", tex_label="\\zeta_{nn}")
    m <= SteadyStateParameter(:ζ_nμ_e,    NaN, description="No description available.", tex_label="\\zeta_{n_{\\mu_e}}")
    m <= SteadyStateParameter(:ζ_nσ_ω,   NaN, description="No description available.", tex_label="\\zeta_{n_{\\sigma_\\omega}}")
end

"""
```
steadystate!(m::Model904)
```

Calculates the model's steady-state values. `steadystate!(m)` must be called whenever
the parameters of `m` are updated.
"""
function steadystate!(m::Model904)
    SIGWSTAR_ZERO = 0.5

    m[:z_star]   = log(1+m[:γ]) + m[:α]/(1-m[:α])*log(m[:Upsilon])
    m[:rstar]    = exp(m[:σ_c]*m[:z_star]) / m[:β]
    m[:Rstarn]   = 100*(m[:rstar]*m[:π_star] - 1)
    m[:r_k_star] = m[:spr]*m[:rstar]*m[:Upsilon] - (1-m[:δ])
    m[:wstar]    = (m[:α]^m[:α] * (1-m[:α])^(1-m[:α]) * m[:r_k_star]^(-m[:α]) / m[:Φ])^(1/(1-m[:α]))
    m[:Lstar]    = 1.
    m[:kstar]    = (m[:α]/(1-m[:α])) * m[:wstar] * m[:Lstar] / m[:r_k_star]
    m[:kbarstar] = m[:kstar] * (1+m[:γ]) * m[:Upsilon]^(1 / (1-m[:α]))
    m[:istar]    = m[:kbarstar] * (1-((1-m[:δ])/((1+m[:γ]) * m[:Upsilon]^(1/(1-m[:α])))))
    m[:ystar]    = (m[:kstar]^m[:α]) * (m[:Lstar]^(1-m[:α])) / m[:Φ]
    m[:cstar]    = (1-m[:g_star])*m[:ystar] - m[:istar]
    m[:wl_c]     = (m[:wstar]*m[:Lstar])/(m[:cstar]*m[:λ_w])

    # FINANCIAL FRICTIONS ADDITIONS
    # solve for σ_ω_star and zω_star
    zω_star = quantile(Normal(), m[:Fω].scaledvalue)

    σ_ω_star = SIGWSTAR_ZERO
    try
        σ_ω_star = fzero(sigma -> ζ_spb_fn(zω_star, sigma, m[:spr]) - m[:ζ_spb], 0.5)
    catch ex
        σ_ω_star = SIGWSTAR_ZERO
        if !isa(ex, ConvergenceFailed)
            rethrow(ex)
        else
            σ_ω_star = SIGWSTAR_ZERO
        end
    end

    # evaluate ωbarstar
    ωbarstar = ω_fn(zω_star, σ_ω_star)

    # evaluate all BGG function elasticities
    Gstar       = G_fn(zω_star, σ_ω_star)
    Γstar       = Γ_fn(zω_star, σ_ω_star)
    dGdω_star   = dG_dω_fn(zω_star, σ_ω_star)
    d2Gdω2star  = d2G_dω2_fn(zω_star, σ_ω_star)
    dΓdω_star   = dΓ_dω_fn(zω_star)
    d2Γdω2star  = d2Γ_dω2_fn(zω_star, σ_ω_star)
    dGdσstar    = dG_dσ_fn(zω_star, σ_ω_star)
    d2Gdωdσstar = d2G_dωdσ_fn(zω_star, σ_ω_star)
    dΓdσstar    = dΓ_dσ_fn(zω_star, σ_ω_star)
    d2Γdωdσstar = d2Γ_dωdσ_fn(zω_star, σ_ω_star)

    # evaluate μ, nk, and Rhostar
    μ_estar     = μ_fn(zω_star, σ_ω_star, m[:spr])
    nkstar      = nk_fn(zω_star, σ_ω_star, m[:spr])
    Rhostar     = 1/nkstar - 1

    # evaluate wekstar and vkstar
    betabar     = exp( (m[:σ_c] -1) * m[:z_star]) / m[:β]
    wekstar     = (1-(m[:γ_star]*betabar))*nkstar - m[:γ_star]*betabar*(m[:spr]*(1-μ_estar*Gstar) - 1)
    vkstar      = (nkstar-wekstar)/m[:γ_star]

    # evaluate nstar and vstar
    m[:nstar]   = nkstar*m[:kbarstar]
    m[:vstar]   = vkstar*m[:kbarstar]

    # a couple of combinations
    ΓμG         = Γstar - μ_estar*Gstar
    ΓμGprime    = dΓdω_star - μ_estar*dGdω_star

    # elasticities wrt ωbar
    ζ_bw       = ζ_bω_fn(zω_star, σ_ω_star, m[:spr])
    ζ_zw       = ζ_zω_fn(zω_star, σ_ω_star, m[:spr])
    ζ_bw_zw    = ζ_bw/ζ_zw

    # elasticities wrt σ_ω
    ζ_bσ_ω      = σ_ω_star * (((1 - μ_estar*dGdσstar/dΓdσstar) /
        (1 - μ_estar*dGdω_star/dΓdω_star) - 1)*dΓdσstar*m[:spr] + μ_estar*nkstar*
            (dGdω_star*d2Γdωdσstar - dΓdω_star*d2Gdωdσstar)/ΓμGprime^2) /
                ((1 - Γstar)*m[:spr] + dΓdω_star/ΓμGprime*(1-nkstar))
    ζ_zσ_ω      = σ_ω_star * (dΓdσstar - μ_estar*dGdσstar) / ΓμG
    m[:ζ_spσ_ω] = (ζ_bw_zw*ζ_zσ_ω - ζ_bσ_ω) / (1-ζ_bw_zw)

    # elasticities wrt μ_e
    ζ_bμ_e      = -μ_estar * (nkstar*dΓdω_star*dGdω_star/ΓμGprime+dΓdω_star*Gstar*m[:spr]) /
        ((1-Γstar)*ΓμGprime*m[:spr] + dΓdω_star*(1-nkstar))
    ζ_zμ_e      = -μ_estar*Gstar/ΓμG
    m[:ζ_spμ_e] = (ζ_bw_zw*ζ_zμ_e - ζ_bμ_e) / (1-ζ_bw_zw)

    # some ratios/elasticities
    Rkstar     = m[:spr]*m[:π_star]*m[:rstar] # (r_k_star+1-δ)/Upsilon*π_star
    ζ_gw       = dGdω_star/Gstar*ωbarstar
    ζ_Gσ_ω     = dGdσstar/Gstar*σ_ω_star

    # elasticities for the net worth evolution
    m[:ζ_nRk]    = m[:γ_star]*Rkstar/m[:π_star]/exp(m[:z_star])*(1+Rhostar)*(1 - μ_estar*Gstar*(1 - ζ_gw/ζ_zw))
    m[:ζ_nR]     = m[:γ_star]*betabar*(1+Rhostar)*(1 - nkstar + μ_estar*Gstar*m[:spr]*ζ_gw/ζ_zw)
    m[:ζ_nqk]    = m[:γ_star]*Rkstar/m[:π_star]/exp(m[:z_star])*(1+Rhostar)*(1 - μ_estar*Gstar*(1+ζ_gw/ζ_zw/Rhostar)) - m[:γ_star]*betabar*(1+Rhostar)
    m[:ζ_nn]     = m[:γ_star]*betabar + m[:γ_star]*Rkstar/m[:π_star]/exp(m[:z_star])*(1+Rhostar)*μ_estar*Gstar*ζ_gw/ζ_zw/Rhostar
    m[:ζ_nμ_e]   = m[:γ_star]*Rkstar/m[:π_star]/exp(m[:z_star])*(1+Rhostar)*μ_estar*Gstar*(1 - ζ_gw*ζ_zμ_e/ζ_zw)
    m[:ζ_nσ_ω]  = m[:γ_star]*Rkstar/m[:π_star]/exp(m[:z_star])*(1+Rhostar)*μ_estar*Gstar*(ζ_Gσ_ω-ζ_gw/ζ_zw*ζ_zσ_ω)

    return m
end

function model_settings!(m::Model904)
    default_settings!(m)

    # Anticipated shocks
    m <= Setting(:n_anticipated_shocks, 0)
    m <= Setting(:n_anticipated_shocks_padding, 0)

    # Conditional data variables
    m <= Setting(:cond_semi_names, [:obs_nominalrate])
    m <= Setting(:cond_full_names, [:obs_gdp, :obs_nominalrate])

    # Forecast
    m <= Setting(:forecast_pseudoobservables, false)
end

"""
```
parameter_groupings(m::Model904)
```

Returns an `OrderedDict{String, Vector{Parameter}}` mapping descriptions of
parameter groupings (e.g. \"Policy Parameters\") to vectors of
`Parameter`s. This dictionary is passed in as a keyword argument to
`prior_table`.
"""
function parameter_groupings(m::Model904)
    policy     = [:ψ1, :ψ2, :ψ3, :ρ, :ρ_rm, :σ_r_m]
    sticky     = [:ζ_p, :ι_p, :ϵ_p, :ζ_w, :ι_w, :ϵ_w]
    other_endo = [:γ, :α, :β, :σ_c, :h, :ν_l, :δ, :Φ, :S′′, :ppsi,
                  :Lmean, :λ_w, :π_star, :g_star]
    financial  = [:Fω, :spr, :ζ_spb, :γ_star]
    processes  = [:ρ_g, :ρ_b, :ρ_μ, :ρ_z, :ρ_σ_w, :ρ_π_star, :ρ_λ_f, :ρ_λ_w, :η_λ_f, :η_λ_w,
                  :σ_g, :σ_b, :σ_μ, :σ_z, :σ_σ_ω, :σ_π_star, :σ_λ_f, :σ_λ_w, :η_gz]

    all_keys     = Vector[policy, sticky, other_endo, financial, processes]
    all_params   = map(keys -> [m[θ]::Parameter for θ in keys], all_keys)
    descriptions = ["Policy", "Nominal Rigidities",
                    "Other Endogenous Propagation and Steady State",
                    "Financial Frictions", "Exogenous Process"]

    groupings = OrderedDict{String, Vector{Parameter}}(zip(descriptions, all_params))

    # Ensure no parameters missing
    incl_params = vcat(collect(values(groupings))...)
    excl_params = [m[θ] for θ in [:Upsilon, :ρ_μ_e, :ρ_γ, :σ_μ_e, :σ_γ]]
    @assert isempty(setdiff(m.parameters, vcat(incl_params, excl_params)))
    @assert isempty(setdiff(vcat(incl_params, excl_params), m.parameters))

    return groupings
end

"""
```
shock_groupings(m::Model904)
```

Returns a `Vector{ShockGroup}`, which must be passed in to
`plot_shock_decomposition`. See `?ShockGroup` for details.
"""
function shock_groupings(m::Model904)
    gov = ShockGroup("g", [:g_sh], RGB(0.70, 0.13, 0.13)) # firebrick
    bet = ShockGroup("b", [:b_sh], RGB(0.3, 0.3, 1.0))
    fin = ShockGroup("FF", [:γ_sh, :μ_e_sh, :σ_ω_sh], RGB(0.29, 0.0, 0.51)) # indigo
    tfp = ShockGroup("z", [:z_sh], RGB(1.0, 0.55, 0.0)) # darkorange
    pmu = ShockGroup("p-mkp", [:λ_f_sh], RGB(0.60, 0.80, 0.20)) # yellowgreen
    wmu = ShockGroup("w-mkp", [:λ_w_sh], RGB(0.0, 0.5, 0.5)) # teal
    pol = ShockGroup("pol", vcat([:rm_sh], [Symbol("rm_shl$i") for i = 1:n_anticipated_shocks(m)]),
                     RGB(1.0, 0.84, 0.0)) # gold
    pis = ShockGroup("pi-LR", [:π_star_sh], RGB(1.0, 0.75, 0.793)) # pink
    mei = ShockGroup("mu", [:μ_sh], :cyan)
    det = ShockGroup("dt", [:dettrend], :gray40)

    return [gov, bet, fin, tfp, pmu, wmu, pol, pis, mei, det]
end
