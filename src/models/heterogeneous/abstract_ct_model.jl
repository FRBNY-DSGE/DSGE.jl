"""
```
AbstractCTModel{T} <: AbstractModel{T}
```

The AbstractCTModel is defined as a subtype of AbstractModel to accomodate the
numerical methods and procedures specific to continuous time models.

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

* `spec::String`: The model specification identifier, \"krusell_smith\", cached
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
abstract type AbstractCTModel{T} <: AbstractModel{T} end
#=
    parameters::ParameterVector{T}                          # vector of all time-invariant model parameters
    steady_state::ParameterVector{T}                        # model steady-state values
    keys::OrderedDict{Symbol,Int}                           # human-readable names for all the model
                                                            # parameters and steady-states

    endogenous_states::OrderedDict{Symbol,Vector{Int}}      # these fields used to create matrices in the
    exogenous_shocks::OrderedDict{Symbol,Int}               # measurement and equilibrium condition equations.
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
=#