"""
```
AnSchorfheide{T} <: AbstractModel{T}
```
The `AnSchorfheide` type defines the structure of the model described in 'Bayesian Estimation of DSGE Models' by Edward P. Herbst and Frank Schorfheide

### Fields

#### Parameters and Steady-States
* `parameters::Vector{AbstractParameter}`: Vector of all time-invariant model parameters.

* `steady_state::Vector{AbstractParameter}`: Model steady-state values, computed as a function of elements of
  `parameters`.

* `keys::Dict{Symbol,Int}`: Maps human-readable names for all model parameters and
  steady-states to their indices in `parameters` and `steady_state`.

#### Inputs to Measurement and Equilibrium Condition Equations

The following fields are dictionaries that map human-readible names to row and column
indices in the matrix representations of of the measurement equation and equilibrium
conditions.

* `endogenous_states::Dict{Symbol,Int}`: Maps each state to a column in the measurement and
  equilibrium condition matrices.

* `exogenous_shocks::Dict{Symbol,Int}`: Maps each shock to a column in the measurement and
  equilibrium condition matrices.

* `expected_shocks::Dict{Symbol,Int}`: Maps each expected shock to a column in the
  measurement and equilibrium condition matrices.

* `equilibrium_conditions::Dict{Symbol,Int}`: Maps each equlibrium condition to a row in the
  model's equilibrium condition matrices.

* `endogenous_states_augmented::Dict{Symbol,Int}`: Maps lagged states to their columns in
  the measurement and equilibrium condition equations. These are added after Gensys solves the
  model.

* `observables::Dict{Symbol,Int}`: Maps each observable to a row in the model's measurement
  equation matrices.

#### Model Specifications and Settings

* `spec::AbstractString`: The model specification identifier, \"AnSchorfheide\", cached here for
  filepath computation.

* `subspec::AbstractString`: The model subspecification number, indicating that some
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

* `data_series::Dict{Symbol,Vector{Symbol}}`: A dictionary that
  stores data sources (keys) and lists of series mnemonics
  (values). DSGE.jl will fetch data from the Federal Reserve Bank of
  St. Louis's FRED database; all other data must be downloaded by the
  user. See `load_data` for further details.
"""
type AnSchorfheide{T} <: AbstractModel{T}
    parameters::ParameterVector{T}                  # vector of all time-invariant model parameters
    steady_state::ParameterVector{T}                # model steady-state values
    keys::Dict{Symbol,Int}                          # human-readable names for all the model
                                                    # parameters and steady-n_states

    endogenous_states::Dict{Symbol,Int}             # these fields used to create matrices in the
    exogenous_shocks::Dict{Symbol,Int}              # measurement and equilibrium condition equations.
    expected_shocks::Dict{Symbol,Int}               #
    equilibrium_conditions::Dict{Symbol,Int}        #
    endogenous_states_augmented::Dict{Symbol,Int}   #
    observables::Dict{Symbol,Int}                   #

    spec::ASCIIString                               # Model specification number (eg "AnSchorfheide")
    subspec::ASCIIString                            # Model subspecification (eg "ss0")
    settings::Dict{Symbol,Setting}                  # Settings/flags for computation
    test_settings::Dict{Symbol,Setting}             # Settings/flags for testing mode
    rng::MersenneTwister                            # Random number generator
    testing::Bool                                   # Whether we are in testing mode or not
    data_series::Dict{Symbol,Vector{Symbol}}       # Keys = data sources, values = vector of series mnemonics
    data_transforms::OrderedDict{Symbol,Function}  # functions to transform raw data into input matrix
end

description(m::AnSchorfheide) = "Julia implementation of model defined in 'Bayesian Estimation of DSGE Models' by Edward P. Herbst and Frank Schorfheide: AnSchorfheide, $(m.subspec)"

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
        :y_t, :π_t, :R_t, :y1_t, :g_t, :z_t, :Ey_t1, :Eπ_t1])

    # Exogenous shocks
    exogenous_shocks = collect([
        :z_sh, :g_sh, :R_sh])

    # Expectations shocks
    expected_shocks = collect([
        :Ey_sh, :Eπ_sh])

    # Equilibrium conditions
    equilibrium_conditions = collect([
        :eq_euler, :nk_pcurve, :mp_rule, :shock_1, :shock_2, :shock_3, :shock_4, :shock_5])

    # Additional states added after solving model
    # Lagged states and observables measurement error
    endogenous_states_augmented = []

    # Measurement equation observables
    observables = collect([
        :obs_gdp,              # quarterly output growth
        :obs_corepce,          # inflation (core PCE)
        :obs_ffr])             # federal funds rate

    for (i,k) in enumerate(endogenous_states);            m.endogenous_states[k]            = i end
    for (i,k) in enumerate(exogenous_shocks);             m.exogenous_shocks[k]             = i end
    for (i,k) in enumerate(expected_shocks);              m.expected_shocks[k]              = i end
    for (i,k) in enumerate(equilibrium_conditions);       m.equilibrium_conditions[k]       = i end
    for (i,k) in enumerate(endogenous_states);            m.endogenous_states[k]            = i end
    for (i,k) in enumerate(endogenous_states_augmented); m.endogenous_states_augmented[k] = i+length(endogenous_states) end
    for (i,k) in enumerate(observables);                  m.observables[k]                  = i end
end


function AnSchorfheide(subspec::AbstractString="ss0")

    # Model-specific specifications
    spec               = split(basename(@__FILE__),'.')[1]
    subspec            = subspec
    settings           = Dict{Symbol,Setting}()
    test_settings      = Dict{Symbol,Setting}()
    rng                = MersenneTwister()
    testing            = false

    fred_series        = [:GDPC1, :PCEPILFE, :DFF]#, :CNP16OV]

    data_series        = Dict(:us =>[:obs_gdp, :obs_corepce, :obs_ffr]) # Dict(:fred => fred_series)
    data_transforms    = OrderedDict{Symbol,Function}()

    # initialize empty model
    m = AnSchorfheide{Float64}(
            # model parameters and steady state values
            Vector{AbstractParameter{Float64}}(), Vector{Float64}(), Dict{Symbol,Int}(),

            # model indices
            Dict{Symbol,Int}(), Dict{Symbol,Int}(), Dict{Symbol,Int}(), Dict{Symbol,Int}(), Dict{Symbol,Int}(), Dict{Symbol,Int}(),

            spec,
            subspec,
            settings,
            test_settings,
            rng,
            testing,
            data_series,
            data_transforms)

    # Set settings
    settings_AnSchorfheide!(m)
    default_test_settings!(m)

    # Set data transformations
    init_data_transforms!(m)

    # Initialize parameters
    m <= parameter(:τ,      1.9937, (1e-20, 1e5), (1e-20, 1e5),   DSGE.Exponential(),     GammaAlt(2., 0.5),         fixed=false,
                   description="The intertemporal elasticity of substitution is 1/τ.",
                   tex_label="\\tau")
    m <= parameter(:κ,      0.7306, (1e-20, 1e1), (1e-20, 1e1),   DSGE.SquareRoot(),     Uniform(0,1),         fixed=false,
                   description="",
                   tex_label="\\kappa")
    m <= parameter(:ψ_1,      1.1434, (1 + 1e-20, 1e5), (1+1e-20, 1e5),   DSGE.Exponential(), GammaAlt(1.5, 0.25),         fixed=false,
                   description="",
                   tex_label="\\psi_1")
    m <= parameter(:ψ_2,      0.4536, (1e-20, 1e5), (1e-20, 1e5),   DSGE.Exponential(),     GammaAlt(0.5, 0.25),         fixed=false,
                   description="",
                   tex_label="\\psi_2")
    m <= parameter(:rA,      0.0313, (1e-20, 1e5), (1e-20, 1e5),   DSGE.Exponential(),     GammaAlt(0.5, 0.5),         fixed=false,
                   description="β (discount factor) =  1/(1+ rA/400)",
                   tex_label="\\rA")
    m <= parameter(:π,      8.1508, (1e-20, 1e5), (1e-20, 1e5),   DSGE.Exponential(),     GammaAlt(7., 2.),         fixed=false,
                   description="Target inflation rate",
                   tex_label="\\pi")
    m <= parameter(:γ_Q,      1.5, (1e-20, 1e5), (1e-20, 1e5),   DSGE.Exponential(),     Normal(0.40, 0.20),         fixed=false,

                   description="Growth rate of technology",
                   tex_label="\\gamma_Q")
    m <= parameter(:ρ_R,      0.3847, (1e-20, 1-1e-7), (1e-20, 1-1e-7),   DSGE.SquareRoot(),     Uniform(0,1),         fixed=false,
                   description="AR(1) coefficient on interest rate",
                   tex_label="\\rho_R")
    m <= parameter(:ρ_g,      0.3777, (1e-20, 1-1e-7), (1e-20, 1-1e-7),   DSGE.SquareRoot(),     Uniform(0,1),         fixed=false,
                   description="AR(1) coefficient on g_t = 1/(1 - ζ_t), where ζ_t is government spending as a fraction of output.",
                   tex_label="\\rho_g")
    m <= parameter(:ρ_z,      0.9579, (1e-20, 1-1e-7), (1e-20, 1-1e-7),   DSGE.SquareRoot(),     Uniform(0,1),         fixed=false,
                   description="",
                   tex_label="\\rho_z")
    m <= parameter(:σ_R,      0.4900, (1e-20, 1e5), (1e-20, 1e5),   DSGE.Exponential(),     DSGE.RootInverseGamma(4, .4),         fixed=false,
                   description="",
                   tex_label="\\sigma_R")
    m <= parameter(:σ_g,      1.4594, (1e-20, 1e5), (1e-20, 1e5),   DSGE.Exponential(),     DSGE.RootInverseGamma(4, 1.),         fixed=false,
                   description="",
                   tex_label="\\sigma_g")
    m <= parameter(:σ_z,      0.9247, (1e-20, 1e5), (1e-20, 1e5),   DSGE.Exponential(),     DSGE.RootInverseGamma(4, 0.5),         fixed=false,
                   description="",
                   tex_label="\\sigma_z")

    init_model_indices!(m)
    init_subspec!(m)
    steadystate!(m)
    return m
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

function settings_AnSchorfheide!(m::AnSchorfheide)
    default_settings!(m)
end
# """
# ```
# init_data_transforms!(m::AnSchorfheide)
# ```

# This function initializes a dictionary of functions that map series read in in levels to the
# appropriate transformed value. At the time that the functions are initialized, data is not
# itself in memory. These functions are model-specific because they assume that certain series
# are available. The keys of data transforms should match exactly the keys of `m.observables`.
# """
# function init_data_transforms!(m::AnSchorfheide)

#     m.data_transforms[:obs_gdp] = function (levels)
#         # FROM: Level of real GDP (from FRED)
#         # TO: Quarter-to-quter percent change of real GDP
#         oneqtrpercentchange(levels[:GDPC1])
#     end

#     m.data_transforms[:obs_corepce] = function (levels)
#         # FROM: Core PCE index (from FRED)
#         # TO: Quarter-to-quter percent change of core PCE, i.e. quarterly inflation
#         oneqtrpercentchange(levels[:PCEPILFE])
#     end

#     m.data_transforms[:obs_ffr] = function (levels)
#         # FROM: Nominal effective federal funds rate (aggregate daily data at a
#         #       quarterly frequency at an annual rate)
#         # TO:   Nominal effective fed funds rate, at a quarterly rate

#         annualtoquarter(levels[:DFF])

#     end
# end


function transform_data(m::AnSchorfheide, levels::DataFrame; verbose::Symbol = :low)
    return levels
end

function init_data_transforms!(m::AnSchorfheide)

    m.data_transforms[:obs_gdp] = function (levels)
    end

    m.data_transforms[:obs_corepce] = function (levels)
    end

    m.data_transforms[:obs_ffr] = function (levels)
    end
end
