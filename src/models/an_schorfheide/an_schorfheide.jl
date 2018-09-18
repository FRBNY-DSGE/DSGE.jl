"""
```
AnSchorfheide{T} <: AbstractModel{T}
```

The `AnSchorfheide` type defines the structure of the simple New Keynesian DSGE
model described in 'Bayesian Estimation of DSGE Models' by Sungbae An and Frank
Schorfheide.

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

* `pseudo_observable_mappings::OrderedDict{Symbol,PseudoObservable}`: A
  dictionary that stores names and transformations to/from model units. See
  `PseudoObservable` for further details.
"""
mutable struct AnSchorfheide{T} <: AbstractModel{T}
    parameters::ParameterVector{T}                         # vector of all time-invariant model parameters
    steady_state::ParameterVector{T}                       # model steady-state values
    keys::OrderedDict{Symbol,Int}                          # human-readable names for all the model
                                                           # parameters and steady-states

    endogenous_states::OrderedDict{Symbol,Int}             # these fields used to create matrices in the
    exogenous_shocks::OrderedDict{Symbol,Int}              # measurement and equilibrium condition equations.
    expected_shocks::OrderedDict{Symbol,Int}               #
    equilibrium_conditions::OrderedDict{Symbol,Int}        #
    endogenous_states_augmented::OrderedDict{Symbol,Int}   #
    observables::OrderedDict{Symbol,Int}                   #
    pseudo_observables::OrderedDict{Symbol,Int}            #

    spec::String                                           # Model specification number (eg "m990")
    subspec::String                                        # Model subspecification (eg "ss0")
    settings::Dict{Symbol,Setting}                         # Settings/flags for computation
    test_settings::Dict{Symbol,Setting}                    # Settings/flags for testing mode
    rng::MersenneTwister                                   # Random number generator
    testing::Bool                                          # Whether we are in testing mode or not

    observable_mappings::OrderedDict{Symbol, Observable}
    pseudo_observable_mappings::OrderedDict{Symbol, PseudoObservable}
end

description(m::AnSchorfheide) = "Julia implementation of model defined in 'Bayesian Estimation of DSGE Models' by Sungbae An and Frank Schorfheide: AnSchorfheide, $(m.subspec)"

"""
`init_model_indices!(m::AnSchorfheide)`

Arguments:
`m:: AnSchorfheide`: a model object

Description:
Initializes indices for all of `m`'s states, shocks, and equilibrium conditions.
"""
function init_model_indices!(m::AnSchorfheide)
    # Endogenous states
    endogenous_states = collect([
        :y_t, :π_t, :R_t, :y_t1, :g_t, :z_t, :Ey_t1, :Eπ_t1])

    # Exogenous shocks
    exogenous_shocks = collect([
        :z_sh, :g_sh, :rm_sh])

    # Expectations shocks
    expected_shocks = collect([
        :Ey_sh, :Eπ_sh])

    # Equilibrium conditions
    equilibrium_conditions = collect([
        :eq_euler, :eq_phillips, :eq_mp, :eq_y_t1, :eq_g, :eq_z, :eq_Ey, :eq_Eπ])

    # Additional states added after solving model
    # Lagged states and observables measurement error
    endogenous_states_augmented = []

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


function AnSchorfheide(subspec::String="ss0";
                       custom_settings::Dict{Symbol, Setting} = Dict{Symbol, Setting}(),
                       testing = false)

    # Model-specific specifications
    spec               = split(basename(@__FILE__),'.')[1]
    subspec            = subspec
    settings           = Dict{Symbol,Setting}()
    test_settings      = Dict{Symbol,Setting}()
    rng                = MersenneTwister(0)

    # initialize empty model
    m = AnSchorfheide{Float64}(
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
    settings_an_schorfheide!(m)
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
init_parameters!(m::AnSchorfheide)
```

Initializes the model's parameters, as well as empty values for the steady-state
parameters (in preparation for `steadystate!(m)` being called to initialize
those).
"""
function init_parameters!(m::AnSchorfheide)
    # Initialize parameters
    m <= parameter(:τ, 1.9937, (1e-20, 1e5), (1e-20, 1e5), DSGE.Exponential(), GammaAlt(2., 0.5), fixed=false,
                   description="τ: The inverse of the intemporal elasticity of substitution.",
                   tex_label="\\tau")

    m <= parameter(:κ, 0.7306, (1e-20, 1e1), (1e-20, 1e1), DSGE.SquareRoot(), Uniform(0,1), fixed=false,
                   description="κ: Composite parameter in New Keynesian Phillips Curve.",
                   tex_label="\\kappa")

    m <= parameter(:ψ_1, 1.1434, (1 + 1e-20, 1e5), (1+1e-20, 1e5), DSGE.Exponential(), GammaAlt(1.5, 0.25), fixed=false,
                   description="ψ_1: The weight on inflation in the monetary policy rule.",
                   tex_label="\\psi_1")
    m <= parameter(:ψ_2, 0.4536, (1e-20, 1e5), (1e-20, 1e5), DSGE.Exponential(), GammaAlt(0.5, 0.25), fixed=false,
                   description="ψ_2: The weight on the output gap in the monetary policy rule.",
                   tex_label="\\psi_2")

    m <= parameter(:rA, 0.0313, (1e-20, 1e5), (1e-20, 1e5), DSGE.Exponential(), GammaAlt(0.5, 0.5), fixed=false,
                   description="rA: β (discount factor) = 1/(1+ rA/400).",
                   tex_label="\\rA")

    m <= parameter(:π_star, 8.1508, (1e-20, 1e5), (1e-20, 1e5), DSGE.Exponential(), GammaAlt(7., 2.), fixed=false,
                   description="π_star: Target inflation rate.",
                   tex_label="\\pi^*")

    m <= parameter(:γ_Q, 1.5, (1e-20, 1e5), (1e-20, 1e5), DSGE.Exponential(), Normal(0.40, 0.20), fixed=false,

                   description="γ_Q: Steady state growth rate of technology.",
                   tex_label="\\gamma_Q")

    m <= parameter(:ρ_R, 0.3847, (1e-20, 1-1e-7), (1e-20, 1-1e-7), DSGE.SquareRoot(), Uniform(0,1), fixed=false,
                   description="ρ_R: AR(1) coefficient on interest rate.",
                   tex_label="\\rho_R")

    m <= parameter(:ρ_g, 0.3777, (1e-20, 1-1e-7), (1e-20, 1-1e-7), DSGE.SquareRoot(), Uniform(0,1), fixed=false,
                   description="ρ_g: AR(1) coefficient on g_t = 1/(1 - ζ_t), where ζ_t is government spending as a fraction of output.",
                   tex_label="\\rho_g")

    m <= parameter(:ρ_z, 0.9579, (1e-20, 1-1e-7), (1e-20, 1-1e-7), DSGE.SquareRoot(), Uniform(0,1), fixed=false,
                   description="ρ_z: AR(1) coefficient on shocks to the technology growth rate.",
                   tex_label="\\rho_z")

    m <= parameter(:σ_R, 0.4900, (1e-20, 1e5), (1e-20, 1e5), DSGE.Exponential(), DSGE.RootInverseGamma(4, .4), fixed=false,
                   description="σ_R: Standard deviation of shocks to the nominal interest rate.",
                   tex_label="\\sigma_R")

    m <= parameter(:σ_g, 1.4594, (1e-20, 1e5), (1e-20, 1e5), DSGE.Exponential(), DSGE.RootInverseGamma(4, 1.), fixed=false,
                   description="σ_g: Standard deviation of shocks to the government spending process.",
                   tex_label="\\sigma_g")

    m <= parameter(:σ_z, 0.9247, (1e-20, 1e5), (1e-20, 1e5), DSGE.Exponential(), DSGE.RootInverseGamma(4, 0.5), fixed=false,
                   description="σ_z: Standard deviation of shocks to the technology growth rate process.",
                   tex_label="\\sigma_z")

    m <= parameter(:e_y, 0.20*0.579923, fixed=true,
                   description="e_y: Measurement error on GDP growth.",
                   tex_label="e_y")

    m <= parameter(:e_π, 0.20*1.470832, fixed=true,
                   description="e_π: Measurement error on inflation.",
                   tex_label="e_\\pi")

    m <= parameter(:e_R, 0.20*2.237937, fixed=true,
                   description="e_R: Measurement error on the interest rate.",
                   tex_label="e_R")
end

"""
```
steadystate!(m::AnSchorfheide)
```

Calculates the model's steady-state values. `steadystate!(m)` must be called whenever
the parameters of `m` are updated.
"""
function steadystate!(m::AnSchorfheide)
    return m
end

function settings_an_schorfheide!(m::AnSchorfheide)
    default_settings!(m)

    # Data
    m <= Setting(:data_id, 0, "Dataset identifier")
    m <= Setting(:cond_full_names, [:obs_gdp, :obs_nominalrate],
        "Observables used in conditional forecasts")
    m <= Setting(:cond_semi_names, [:obs_nominalrate],
        "Observables used in semiconditional forecasts")

    # Metropolis-Hastings
    m <= Setting(:mh_cc, 0.27,
                 "Jump size for Metropolis-Hastings (after initialization)")

    # Forecast
    m <= Setting(:use_population_forecast, true,
                 "Whether to use population forecasts as data")
    m <= Setting(:forecast_zlb_value, 0.13,
        "Value of the zero lower bound in forecast periods, if we choose to enforce it")
end

function shock_groupings(m::AnSchorfheide)
    gov = ShockGroup("g", [:g_sh], RGB(0.70, 0.13, 0.13)) # firebrick
    tfp = ShockGroup("z", [:z_sh], RGB(1.0, 0.55, 0.0)) # darkorange
    pol = ShockGroup("pol", vcat([:rm_sh], [Symbol("rm_shl$i") for i = 1:n_anticipated_shocks(m)]),
                     RGB(1.0, 0.84, 0.0)) # gold
    det = ShockGroup("dt", [:dettrend], :gray40)

    return [gov, tfp, pol, det]
end
