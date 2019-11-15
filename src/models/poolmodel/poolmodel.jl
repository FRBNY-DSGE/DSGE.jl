"""
```
PoolModel{T}
```

Implements several model averaging methods to pool
the predictions of different models. Currently, `PoolModel`
only works for pooling two different models, although
it can be extended to more models. Confer with the paper for details.

The available methods are

* dynamic: weights on pooled models evolve over time and are
           treated as a stochastic process.

* static:  weights on pooled models are assumed time-invariant.

* equal:   weights are set to 1/2.

* bma:     weights are updated according to Bayesian model averaging.

The default method is dynamic, which implements the Dynamic Pools method
developed by Del Negro et al. (2016). To choose another method,
use the keyword `weight_type::Symbol`, e.g.

```
pm = PoolModel(weight_type = :static) # creates a static PoolModel
```

The `weight_type` is stored as a setting, so users can retrieve
at any point by using `get_setting`.

### Fields

#### Parameters
* `parameters::Vector{AbstractParameter}`: Vector of all time-invariant hyperparameters
    for the chosen method of model averaging.

* `keys::OrderedDict{Symbol,Int}`: Maps human-readable names for predictive densities
    of pooled models.

#### Inputs to Measurement and Equilibrium Condition Equations

* `model::OrderedDict{Symbol,AbstractDSGEModel}`: Maps name to its underlying model
  object.

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
mutable struct PoolModel{T} <: AbstractDSGEModel{T}
    parameters::ParameterVector{T}                         # vector of all time-invariant model parameters
    keys::OrderedDict{Symbol,Int}                          # human-readable names for all the model
    observables::OrderedDict{Symbol,Int}                   # Model names to observables (predictive densities)

    spec::String                                           # Model specification number (eg "m990")
    subspec::String                                        # Model subspecification (eg "ss0")
    settings::Dict{Symbol,Setting}                         # Settings/flags for computation
    test_settings::Dict{Symbol,Setting}                    # Settings/flags for testing mode
    rng::MersenneTwister                                   # Random number generator
    testing::Bool                                          # Whether we are in testing mode or not

    observable_mappings::OrderedDict{Symbol, Observable}
end

description(m::PoolModel) = "Julia implementation of prediction pool methods defined in 'Dynamic prediction pools: An investigation of financial frictions and forecasting performance' by Marco Del Negro, Raiden B. Hasegawa, and Frank Schorfheide: PoolModel, $(m.subspec)"

function init_model_indices!(m::PoolModel)
    # Observables
    observables = keys(m.observable_mappings)

    # Make indices
    for (i,k) in enumerate(observables);       m.observables[k]       = i end
end

function PoolModel(subspec::String="ss2";
                   custom_settings::Dict{Symbol,Setting} = Dict{Symbol,Setting}(),
                   testing = false, verbose::Symbol = :low,
                   weight_type::Symbol = :dynamic) where T<:AbstractFloat

    # Model-specific specifications
    spec               = split(basename(@__FILE__),'.')[1]
    subspec            = subspec
    settings           = Dict{Symbol,Setting}()
    test_settings      = Dict{Symbol,Setting}()
    rng                = MersenneTwister(0)

    # initialize empty model
    m = PoolModel{Float64}(
        # model parameters
        Vector{AbstractParameter{Float64}}(), OrderedDict{Symbol,Int}(),

        # Observable indices
        OrderedDict{Symbol,Int}(),

        # settings
        spec,
        subspec,
        settings,
        test_settings,
        rng,
        testing,
        OrderedDict{Symbol,Observable}())

    # Set settings
    model_settings!(m; weight_type = weight_type)
    default_test_settings!(m)
    for custom_setting in values(custom_settings)
        m <= custom_setting
    end

    # Set observable and pseudo-observable transformations
    init_observable_mappings!(m)

    # Initialize parameters
    init_parameters!(m)

    # Initialize model indices and subspec
    init_model_indices!(m)
    init_subspec!(m)

    return m
end

"""
```
init_parameters!(m::PoolModel)
```

Initializes the model's parameters, as well as empty values for the steady-state
parameters (in preparation for `steadystate!(m)` being called to initialize
those).
"""
function init_parameters!(m::PoolModel)
    # Initialize parameters
    weight_type = get_setting(m, :weight_type)
    if weight_type == :dynamic
        m <= parameter(:ρ, 0.8, (1e-5,0.999), (1e-5,0.999), SquareRoot(), Uniform(0.,1.),
                       fixed = false,
                       description="ρ: persistence of AR processing underlying λ.",
                       tex_label="\\rho")
        m <= parameter(:μ, 0., fixed = true,
                       description="μ: drift of AR processing underlying λ.",
                       tex_label="\\mu")
        m <= parameter(:σ, 1., fixed = true,
                       description="σ: volatility of AR processing underlying λ.",
                       tex_label="\\sigma")
    elseif weight_type == :equal
        m <= parameter(:λ, 0.5, (1e-5,0.999), (1e-5,0.999), SquareRoot(), Uniform(0.,1.),
                       fixed = true, description="λ: weight on model 1's predictive density",
                       tex_label="\\lambda")
    elseif weight_type == :static
        m <= parameter(:λ, 0.5, (1e-5,0.999), (1e-5,0.999), SquareRoot(), Uniform(0.,1.),
                       fixed = false, description="λ: weight on model 1's predictive density",
                       tex_label="\\lambda")
    elseif weight_type == :bma
        m <= parameter(:λ, 0.5, (1e-5,0.999), (1e-5,0.999), SquareRoot(), Uniform(0.,1.),
                       fixed = false, description="λ: weight on model 1's predictive density",
                       tex_label="\\lambda")
    end
end

function model_settings!(m::PoolModel; weight_type::Symbol = :dynamic_weight)
    default_settings!(m)

    # Weight type: dynamic, equal, or static
    if !(weight_type in [:dynamic, :equal, :static, :bma])
        error("Weight type of a PoolModel object must be :dynamic, :equal, :static, or :bma.")
    else
        m <= Setting(:weight_type, weight_type, "How to weight predictive densities")
    end

    # Data
    m <= Setting(:population_mnemonic, Nullable())
    m <= Setting(:data_id, 1922016, "Dataset identifier")
    m <= Setting(:date_presample_start, quartertodate("1992-Q1"))
    m <= Setting(:date_mainsample_start, quartertodate("1992-Q1"))
    m <= Setting(:date_forecast_start, quartertodate("2011-Q3"))

    # SMC estimation
    m <= Setting(:sampling_method, :SMC)
    m <= Setting(:n_particles, 2000)

    # Forecast
    m <= Setting(:use_population_forecast, false,
                 "Whether to use population forecasts as data")

    # Tempered particle filter
    m <= Setting(:fixed_sched, [1.],
                 "schedule for tempering in tpf; leave empty if want adaptive tempering")
    tuning = Dict(:r_star => 2., :c_init => 0.3, :target_accept_rate => 0.4,
              :resampling_method => :systematic, :n_mh_steps => 1,
              :n_particles => 1000, :n_presample_periods => 0,
              :allout => true)
    m <= Setting(:tuning, tuning, "tuning parameters for TPF")

    # MH estimation
    m <= Setting(:mh_cc, 0.15)
    m <= Setting(:mh_cc0, 0.15)
    m <= Setting(:calculate_hessian, false)
    m <= Setting(:reoptimize, false)
    m <= Setting(:n_mh_simulations, 1000)
    m <= Setting(:n_mh_blocks, 1)
    m <= Setting(:mh_thin, 1)
    m <= Setting(:n_mh_burn, 0)
end

#########################################################
# Overloading various functions for PoolModel type
#########################################################
function Base.show(io::IO, m::PoolModel)
    model_str = ""
    n_obs = n_observables(m)
    for i in 1:n_obs
        if i < n_obs
            model_str *= string(m.observable_mappings[get_key(m, :obs, i)].key) * ", "
        else
            model_str *= string(m.observable_mappings[get_key(m, :obs, i)].key) * "\n"
        end
    end
    @printf io "Dynamic Prediction Pool Method\n"
    @printf io "models: %s\n"                 model_str
    @printf io "data vintage:           %s\n" data_vintage(m)
    @printf io "description:\n %s\n"          description(m)
end

"""
```
steadystate!(m::PoolModel)
```

Computes the steady state for a `PoolModel`. Since
the current implementation does not have a useful
reason to compute the steady state, this function does nothing.

### Arguments:
- `m`: the model object
"""
function steadystate!(m::PoolModel)
end
