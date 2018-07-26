import DataStructures: OrderedDict
import DSGE: description, AbstractModel, default_settings!,
             AbstractParameter, ParameterVector, parameter, SteadyStateParameter,
             SquareRoot, BetaAlt, Exponential, RootInverseGamma,
             Setting, get_setting,
             Observable, PseudoObservable,
             init_subspec!,
             DSGE_DATE_FORMAT

"""
```
KrusellSmith{T} <: AbstractModel{T}
```

The `KrusellSmith` type defines the structure of the simple New Keynesian DSGE
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
type KrusellSmith{T} <: AbstractModel{T}
    parameters::ParameterVector{T}                         # vector of all time-invariant model parameters
    steady_state::ParameterVector{T}                       # model steady-state values
    grids::OrderedDict{Symbol,Grid}                        # quadrature grids
    keys::OrderedDict{Symbol,Int}                          # human-readable names for all the model
                                                           # parameters and steady-states

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

description(m::KrusellSmith) = "KrusellSmith, $(m.subspec)"

"""
`init_model_indices!(m::KrusellSmith)`

Arguments:
`m:: KrusellSmith`: a model object

Description:
Initializes indices for all of `m`'s states, shocks, and equilibrium conditions.
"""
function init_model_indices!(m::KrusellSmith)
    # Endogenous states
    endogenous_states = collect([
    # These states corresp. to the following in the original notation
    #    LMP,    LELLP,  KKP,   ZP,    MP,    ELLP,
        :μ′_t1, :l′_t1, :K′_t, :z′_t, :μ′_t, :l′_t])

    # Exogenous shocks
    exogenous_shocks = collect([:z_sh])

    # Equilibrium conditions
    equilibrium_conditions = collect([
        :eq_euler, :eq_kolmogorov_fwd, :eq_μ, :eq_l, :eq_K_law_of_motion, :eq_TFP])

    # Observables
    observables = keys(m.observable_mappings)

    ########################################################################################
    # Setting indices of endogenous_states and equilibrium conditions manually for now
    nw = get_setting(m, :nw)
    endo = m.endogenous_states_unnormalized
    eqconds = m.equilibrium_conditions

    # State variables
    endo[:μ′_t1] = 1:nw
    endo[:l′_t1] = nw+1:2*nw
    endo[:K′_t]  = 2*nw+1:2*nw+1
    endo[:z′_t]  = 2*nw+2:2*nw+2

    # Jump variables
    endo[:μ′_t]  = 2*nw+3:3*nw+2
    endo[:l′_t]  = 3*nw+3:4*nw+2

    eqconds[:eq_euler]              = 1:nw
    eqconds[:eq_kolmogorov_fwd]     = nw+1:2*nw
    eqconds[:eq_μ]                  = 2*nw+1:3*nw
    eqconds[:eq_l]                  = 3*nw+1:4*nw
    eqconds[:eq_K_law_of_motion]    = 4*nw+1:4*nw+1
    eqconds[:eq_TFP]                = 4*nw+2:4*nw+2
    ########################################################################################

    m.endogenous_states = deepcopy(endo)
    for (i,k) in enumerate(exogenous_shocks);            m.exogenous_shocks[k]            = i end
    for (i,k) in enumerate(observables);                 m.observables[k]                 = i end
end

function KrusellSmith(subspec::String="ss0";
                      custom_settings::Dict{Symbol, Setting} = Dict{Symbol, Setting}(),
                      testing = false)

    # Model-specific specifications
    spec               = "KrusellSmith"
    subspec            = subspec
    settings           = Dict{Symbol,Setting}()
    test_settings      = Dict{Symbol,Setting}()
    rng                = MersenneTwister(0)

    # initialize empty model
    m = KrusellSmith{Float64}(
            # model parameters and steady state values
            Vector{AbstractParameter{Float64}}(), Vector{Float64}(),
            # grids and keys
            OrderedDict{Symbol,Grid}(), OrderedDict{Symbol,Int}(),

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

    # Set observable transformations
    init_observable_mappings!(m)

    # Initialize parameters
    init_parameters!(m)

    # Initialize grids
    init_grids!(m)

    init_model_indices!(m)

    # So that the indices of m.endogenous_states reflect the normalization
    normalize_state_indices!(m)

    return m
end

"""
```
init_parameters!(m::KrusellSmith)
```

Initializes the model's parameters, as well as empty values for the steady-state
parameters (in preparation for `steadystate!(m)` being called to initialize
those).
"""
function init_parameters!(m::KrusellSmith)
    # Initialize parameters
    m <= parameter(:β, 0.95, fixed = true,
                   description = "β: Discount rate.", tex_label = "\\beta")
    m <= parameter(:γ, 3.0, fixed = true,
                   description = "γ: CRRA Parameter.", tex_label = "\\gamma")
    m <= parameter(:α, 1.0/3.0, fixed = true,
                   description = "α: Capital Share.", tex_label = "\\alpha")
    m <= parameter(:δ, 0.2, fixed = true,
                   description = "δ: Depreciation.", tex_label = "\\delta")
    m <= parameter(:ρ_z, 0.95, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), BetaAlt(0.5, 0.2), fixed=false,
                   description="ρ_z: AR(1) coefficient in the technology process.",
                   tex_label="\\rho_z")
    m <= parameter(:σ_z, sqrt(.007), (1e-8, 5.), (1e-8, 5.), Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_z: The standard deviation of the process describing the stationary component of productivity.",
                   tex_label="\\sigma_{z}")
    m <= parameter(:μ_s, 0., fixed = true, description = "μ_s: Mu of log normal in income")
    m <= parameter(:σ_s, 0.5, fixed = true,
                   description = "σ_s: Sigma of log normal in income")

    m <= parameter(:e_y, 0.1, fixed = true, description = "e_y: Measurement error on GDP",
                   tex_label = "e_y")

    # Setting steady-state parameters
    nw = get_setting(m, :nw)

    m <= SteadyStateParameterGrid(:lstar, fill(NaN, nw), description = "Steady-state expected discounted
                                  marginal utility of consumption", tex_label = "l_*")
    m <= SteadyStateParameterGrid(:cstar, fill(NaN, nw), description = "Steady-state consumption",
                                  tex_label = "c_*")
    m <= SteadyStateParameterGrid(:μstar, fill(NaN, nw), description = "Steady-state cross-sectional
                                  density of cash on hand", tex_label = "\\mu_*")
    m <= SteadyStateParameter(:Kstar, NaN, description = "Steady-state capital", tex_label = "K_*")
    m <= SteadyStateParameter(:Lstar, NaN, description = "Steady-state labor", tex_label = "L_*")
    m <= SteadyStateParameter(:Wstar, NaN, description = "Steady-state wages", tex_label =
                              "W_*")
    m <= SteadyStateParameter(:Rstar, NaN, description = "Steady-state net-return on
                              capital", tex_label = "R_*")
    m <= SteadyStateParameterGrid(:KFstar, fill(NaN, (nw, nw)), description = "Steady-state Kolmogorov
                                  Forward Equation", tex_label = "KF_*")
end

"""
```
init_grids!(m::KrusellSmith)
```
"""
function init_grids!(m::KrusellSmith)
    wscale  = get_setting(m, :wscale)
    wlo     = get_setting(m, :wlo)
    whi     = get_setting(m, :whi)
    nw      = get_setting(m, :nw)

    slo     = get_setting(m, :slo)
    shi     = get_setting(m, :shi)
    ns      = get_setting(m, :ns)

    grids = OrderedDict()

    grids[:wgrid] = Grid(uniform_quadrature(wscale), wlo, whi, nw, scale = wscale)
    grids[:sgrid] = Grid(curtis_clenshaw_quadrature(2), slo, shi, ns)

    m.grids = grids
end

# """
# ```
# steadystate!(m::KrusellSmith)
# ```

# Calculates the model's steady-state values. `steadystate!(m)` must be called whenever
# the parameters of `m` are updated.
# """
# function steadystate!(m::KrusellSmith)
    # return m
# end

function model_settings!(m::KrusellSmith)
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

    # Anticipated shocks
    m <= Setting(:n_anticipated_shocks, 0,
                 "Number of anticipated policy shocks")
    # m <= Setting(:n_anticipated_shocks_padding, 20,
                 # "Padding for anticipated policy shocks")

    # Number of states and jumps
    # May want to generalize this functionality for other models with multiple
    # distributions
    m <= Setting(:normalize_distr_variables, true, "Whether or not to perform the
                 normalization of the μ distribution in the Klein solution step")

    m <= Setting(:state_indices, 1:4, "Which indices of m.endogenous_states correspond to state
                 variables")
    m <= Setting(:jump_indices, 5:6, "Which indices of m.endogenous_states correspond to jump
                 variables")
    m <= Setting(:n_states, 162 - get_setting(m, :normalize_distr_variables),
                 "Number of state variables, in the true sense (fully
                 backward looking) accounting for the discretization across the grid")
    m <= Setting(:n_jumps, 160 - get_setting(m, :normalize_distr_variables),
                 "Number of jump variables (forward looking) accounting for
                the discretization across the grid")
    m <= Setting(:n_model_states, get_setting(m, :n_states) + get_setting(m, :n_jumps),
                 "Number of 'states' in the state space model. Because backward and forward
                 looking variables need to be explicitly tracked for the Klein solution
                 method, we have n_states and n_jumps")

    # Mollifier setting parameters
    m <= Setting(:In, 0.443993816237631, "Normalizing constant for the mollifier")
    m <= Setting(:smoother, 1.0, "Scale of the molifier (how peaked the mollifier is)
                 smoother = 0 you get uniform dist")
    m <= Setting(:zlo, 0.0, "Second income shock to mollify actual income")
    m <= Setting(:zhi, 2.0, "Upper bound on this process")

    # w: Cash on Hand Grid Setup
    m <= Setting(:nw, 80, "Cash on hand distribution grid points")
    m <= Setting(:wlo, 0.5, "Lower bound on cash on hand")
    m <= Setting(:whi, 15.0, "Upper Bound on cash on hand")
    m <= Setting(:wscale, get_setting(m, :whi) - get_setting(m, :wlo), "Size of the wgrid")

    # s: Skill Distribution/ "Units of effective labor" Grid Setup
    m <= Setting(:ns, 50, "Skill distribution grid points")
    m <= Setting(:slo, 0.5, "Lower bound on skill grid")
    m <= Setting(:shi, 1.5, "Upper bound on skill")
    m <= Setting(:sscale, (get_setting(m, :shi) - get_setting(m, :slo))/2.0, "Size of our s
                 grid, relative to size of [-1, 1]")
end

# For normalizing the states
function normalize_state_indices!(m::AbstractModel)
    endo = m.endogenous_states
    state_indices = get_setting(m, :state_indices)
    jump_indices  = get_setting(m, :jump_indices)
    normalization_factor = get_setting(m, :normalize_distr_variables)

    model_state_keys    = endo.keys
    jump_keys           = endo.keys[jump_indices]

    # This structure assumes that states are always ordered before jumps
    # And that both states and jumps have the same number of variables
    # to be normalized, in this case μ′_t1 and μ′_t
    normalize_states!(endo, normalization_factor, model_state_keys)
    normalize_jumps!(endo, normalization_factor, jump_keys)

    # TO DO: Include assertions to ensure that the indices are all
    # consecutive and that the right factor was subtracted from each distribution object
end

# Shift a UnitRange type down by an increment
# If this UnitRange is the first_range in the group being shifted, then do not
# subtract the increment from the start
# e.g.
# If you have 1:80 as your first range, then shift(1:80, -1; first_range = true)
# should return 1:79
# However, if you have 81:160 as range, that is not your first range, then the
# function should return 80:159
function shift(inds::UnitRange, increment::Int64; first_range::Bool = false)
    if first_range
        return UnitRange(inds.start, inds.stop + increment)
    else
        return UnitRange(inds.start + increment, inds.stop + increment)
    end
end

function normalize_states!(endo::OrderedDict, normalization_factor::Bool,
                           model_state_keys::Vector{Symbol})
    for (i, model_state) in enumerate(model_state_keys)
        inds = endo[model_state]
        if i == 1
            endo[model_state] = shift(inds, -normalization_factor, first_range = true)
        else
            endo[model_state] = shift(inds, -normalization_factor)
        end
    end
end

function normalize_jumps!(endo::OrderedDict, normalization_factor::Bool,
                          jump_keys::Vector{Symbol})
    for (i, jump) in enumerate(jump_keys)
        inds = endo[jump]
        if i == 1
            endo[jump] = shift(inds, -normalization_factor, first_range = true)
        else
            endo[jump] = shift(inds, -normalization_factor)
        end
    end
end
