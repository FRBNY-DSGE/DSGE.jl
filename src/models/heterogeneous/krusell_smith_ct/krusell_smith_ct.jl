"""
```
KrusellSmithCT{T} <: AbstractCTModel{T}
```

The `KrusellSmithCT` mutable struct defines the structure of the Krusell Smith (1998) model
in 'Income and Wealth Heterogeneity in the Macroeconomy'
originally written in MATLAB by Ben Moll and SeHyoun Ahn at
https://sehyoun.com/EXAMPLE_PHACT_KS.html

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
mutable struct KrusellSmithCT{T} <: AbstractCTModel{T}
    parameters::ParameterVector{T}                          # vector of all time-invariant model parameters
    steady_state::ParameterVector{T}                        # model steady-state values
    keys::OrderedDict{Symbol,Int}                           # human-readable names for all the model
                                                            # parameters and steady-states

    endogenous_states::OrderedDict{Symbol,Vector{Int}}      # these fields are used to create matrices in the
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

description(m::KrusellSmithCT) = "Julia implementation of the Krusell Smith (1998) model defined in 'Income and Wealth Heterogeneity in the Macroeconomy' and originally written in MATLAB by Ben Moll and SeHyoun Ahn at https://sehyoun.com/EXAMPLE_PHACT_KS.html: KrusellSmithCT, $(m.subspec)"

"""
`init_model_indices!(m::KrusellSmithCT)`

Arguments:
`m:: KrusellSmithCT`: a model object

Description:
Initializes indices for all of `m`'s states, shocks, and equilibrium conditions.
"""
function init_model_indices!(m::KrusellSmithCT)
    # Endogenous states
    endogenous_states = collect([:value_function, :distribution, :log_aggregate_tfp, :K, :r, :w, :output, :C, :investment])

    # Exogenous shocks
    exogenous_shocks = collect([:aggregate_tfp])

    # Expectations shocks
    expected_shocks = collect([:E_V])

    # Equilibrium conditions
    equilibrium_conditions = collect([:value_function, :distribution, :log_aggregate_tfp, :K, :r, :w, :output, :C, :investment])

    # Additional states added after solving model
    # Lagged states and observables measurement error
    endogenous_states_augmented = []

    # Observables
    observables = keys(m.observable_mappings)

    # Pseudo-observables
    pseudo_observables = keys(m.pseudo_observable_mappings)

    # Assign dimensions
    m.endogenous_states[:value_function]    = collect(1:200)     # size of state space grid
    m.endogenous_states[:distribution]      = collect(201:399)   # size of grid minus one for ...
    m.endogenous_states[:log_aggregate_tfp] = collect([400])     # number of shocks to distribution
    m.endogenous_states[:K]                 = [401]
    m.endogenous_states[:r]                 = [402]
    m.endogenous_states[:w]                 = [403]
    m.endogenous_states[:output]            = [404]
    m.endogenous_states[:C]                 = [405]
    m.endogenous_states[:investment]        = [406]
    m.expected_shocks[:E_V]                 = m.endogenous_states[:value_function]

    for (i,k) in enumerate(observables);             m.observables[k]             = i end
    for (i,k) in enumerate(exogenous_shocks);        m.exogenous_shocks[k]        = i end
    for (i,k) in enumerate(pseudo_observables);      m.pseudo_observables[k]      = i end
end


function KrusellSmithCT(subspec::String="ss0";
                       custom_settings::Dict{Symbol, Setting} = Dict{Symbol, Setting}(),
                       testing = false)

    # Model-specific specifications
    spec               = split(basename(@__FILE__),'.')[1]
    subspec            = subspec
    settings           = Dict{Symbol,Setting}()
    test_settings      = Dict{Symbol,Setting}()
    rng                = MersenneTwister(0)

    # initialize empty model
    m = KrusellSmithCT{Float64}(
            # model parameters and steady state values
            Vector{AbstractParameter{Float64}}(), Vector{Float64}(),
            OrderedDict{Symbol,Int}(),

            # model indices
            OrderedDict{Symbol,Vector{Int}}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Vector{Int}}(), OrderedDict{Symbol,Vector{Int}}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(),

            spec,
            subspec,
            settings,
            test_settings,
            rng,
            testing,
            OrderedDict{Symbol,Observable}(),
            OrderedDict{Symbol,PseudoObservable}())

    # Set parameters
    init_parameters!(m)
    init_model_indices!(m)

    # Set settings
    model_settings!(m)
    default_test_settings!(m)

    for custom_setting in values(custom_settings)
        m <= custom_setting
    end

    # Test to see if needed values are defined, throw error otherwise
    try
        @assert typeof(get_setting(m, :reduce_state_vars)) == Bool
        @assert typeof(get_setting(m, :reduce_v)) == Bool
    catch
        error("Need to specify in settings whether to perform distribution and/or value function reduction")
    end

    if get_setting(m, :reduce_state_vars) || get_setting(m, :reduce_v)
        try
            @assert get_setting(m, :n_jump_vars) > 0
            @assert get_setting(m, :n_state_vars) > 0
            @assert get_setting(m, :n_state_vars_unreduce) >= 0
        catch
            error("Need to enter values in settings for n_jump_vars, n_state_vars, and/or n_state_vars_unreduce in Settings. These indicate the number of jump and aggregate state variables.")
        end
    end

    if get_setting(m, :reduce_state_vars)
        try
            @assert get_setting(m, :krylov_dim) > 0
        catch
            error("Need to specify in settings desired dimension of Krylov subspace, krylov_dim.")
        end
    end

    if get_setting(m, :reduce_v)
        try
            @assert get_setting(m, :n_prior) > 0 && typeof(get_setting(m, :n_prior)) == Int
            @assert get_setting(m, :n_post) > 0 && typeof(get_setting(m, :n_post)) == Int
            @assert in(:knots_dict, keys(m.settings))
            @assert in(:spline_grid, keys(m.settings))
        catch
            error("Need to specify in settings n_prior, n_post, knots_dict and/or spline_grid for value function reduction. See reduction.jl for their uses.")
        end
    end

    # Solve for steady state
    steadystate!(m)
    return m
end

"""
```
init_parameters!(m::KrusellSmithCT)
```

Initializes the model's parameters, as well as empty values for the steady-state
parameters (in preparation for `steadystate!(m)` being called to initialize
those).
"""
function init_parameters!(m::KrusellSmithCT)
    # Initialize parameters
    m <= parameter(:γ, 2., fixed=true,
                   description="γ: The coefficient of relative risk aversion.",
                   tex_label="\\gamma")

    m <= parameter(:ρ, 0.01, fixed=true,
                   description="ρ: The rate of time preference.",
                   tex_label="\\rho")

    m <= parameter(:δ, .025, fixed=true,
                   description="δ: Capital depreciation.",
                   tex_label="\\delta")

    m <= parameter(:α, 1/3, fixed=true,
                   description="α: Capital share.",
                   tex_label="\\alpha")

    m <= parameter(:σ_tfp, .007, (1e-20, 1e5), (1e-20, 1e5), ModelConstructors.Exponential(), RootInverseGamma(4, 0.5),
                   fixed=false, description="σ_tfp: The standard deviation of the TFP shock.",
                   tex_label="\\sigma_tfp")

    m <= parameter(:ρ_tfp, .95, fixed=true,
                   description="ρ_tfp: The quarterly autocorrelation of TFP shock.",
                   tex_label="\\rho_tfp")

    m <= parameter(:λ1, 1/2, fixed=true,
                   description="λ1: The expected duration of unemployment is 2 quarters.",
                   tex_label="\\lambda_1")

    m <= parameter(:target_ur, .07, fixed=true,
                   description="target_ur: Target unemployment rate.",
                   tex_label="target_ur")

    m <= parameter(:target_er, .93, fixed=true,
                   description="target_er: Target employment rate.",
                   tex_label="target_er")

    m <= parameter(:λ2, (m[:λ1] / (1. * m[:target_er] - 0.)) * (1. - 1. * m[:target_er]), fixed=true,
                   description="λ2: The unemployment rate is 7.",
                   tex_label="\\lambda_2")

    m <= parameter(:μ, .15, fixed=true,
                   description="μ: The UI replacement rate is 15.",
                   tex_label="\\mu")

    m <= parameter(:τ, (m[:μ] / 1.) * (m[:λ2] / m[:λ1]) , fixed=true,
                   description="τ: The labor income tax.",
                   tex_label="\\tau")

    # Steady state parameters
    m <= SteadyStateParameterGrid(:V_ss, Vector{Float64}(undef, 0); description =
                                  "Stacked steady-state value function")
    m <= SteadyStateParameterGrid(:gg_ss, Vector{Float64}(undef, 0); description =
                                  "Stacked steady-state distribution across state space")
    m <= SteadyStateParameter(:log_aggregate_tfp_ss, NaN, description =
                              "Steady-state level of aggregate shocks")
    m <= SteadyStateParameter(:K_ss, NaN, description = "Steady-state capital level")
    m <= SteadyStateParameter(:r_ss, NaN, description = "Steady-state real interest rate")
    m <= SteadyStateParameter(:w_ss, NaN, description = "Steady-state wage")
    m <= SteadyStateParameter(:Y_ss, NaN, description = "Steady-state output level")
    m <= SteadyStateParameter(:C_ss, NaN, description = "Steady-state consumption level")
    m <= SteadyStateParameter(:I_ss, NaN, description = "Steady-state investment level")
    m <= SteadyStateParameterGrid(:If_ss, Matrix{Float64}(undef, 0, 0); description =
                                   "Steady-state forward difference upwind scheme")
    m <= SteadyStateParameterGrid(:Ib_ss, Matrix{Float64}(undef, 0, 0); description =
                                  "Steady-state backward difference upwind scheme")
    m <= SteadyStateParameterGrid(:I0_ss, Matrix{Float64}(undef, 0, 0); description =
                                  "Steady-state no drift upwind scheme")
end

 """
 ```
 steadystate!(m::KrusellSmithCT)
 ```

 Calculates the model's steady-state values. `steadystate!(m)` must be called whenever
 the parameters of `m` are updated.
 """
function steadystate!(m::KrusellSmithCT)

    # Read in parameters
    γ::Float64     = m[:γ].value
    ρ::Float64     = m[:ρ].value
    δ::Float64     = m[:δ].value
    α::Float64     = m[:α].value
    σ_tfp::Float64 = m[:σ_tfp].value
    ρ_tfp::Float64 = m[:ρ_tfp].value
    μ::Float64     = m[:μ].value # UI replacement rate
    τ::Float64     = m[:τ].value # labor income tax

    # Read in grids
    I::Int64             = get_setting(m, :I)
    J::Int64             = get_setting(m, :J)
    z::Vector{Float64}   = get_setting(m, :z)
    a::Vector{Float64}   = get_setting(m, :a)
    da::Float64          = get_setting(m, :da)
    dz::Float64          = get_setting(m, :dz)
    amax::Float64        = get_setting(m, :amax)
    amin::Float64        = get_setting(m, :amin)
    z_avg::Float64       = (m[:λ1].value * z[2] + m[:λ2].value * z[1]) / (m[:λ1].value + m[:λ2].value)
    aa::Array{Float64,2} = get_setting(m, :aa)
    zz::Array{Float64,2} = get_setting(m, :zz)

    # A_switch * v = lambda_z * (v(a, z') - v(a, z)) in Eqn 8 of Why Inequality Matters
    # basically captures the expected change in value function due to jumps in z
    A_switch = [-SparseMatrixCSC{Float64}(LinearAlgebra.I, I, I) * m[:λ1].value SparseMatrixCSC{Float64}(LinearAlgebra.I, I, I) * m[:λ1].value; SparseMatrixCSC{Float64}(LinearAlgebra.I, I, I) * m[:λ2].value -SparseMatrixCSC{Float64}(LinearAlgebra.I, I, I) * m[:λ2].value]

    # Read in initial rates
    r0::Float64   = get_setting(m, :r0)
    rmin::Float64 = get_setting(m, :rmin)
    rmax::Float64 = get_setting(m, :rmax)

    # Read in  approximation parameters
    Ir::Int64          = get_setting(m, :Ir)
    Δ::Float64         = get_setting(m, :Δ_HJB)
    crit_S::Float64    = get_setting(m, :crit_S)
    crit_HJB::Float64  = get_setting(m, :crit_HJB)
    maxit_HJB::Float64 = get_setting(m, :maxit_HJB)

    # Initialize guesses and difference matrices
    dVf::Array{Float64,2}       = zeros(Float64, I, 2) # forward value function difference matrix
    dVb::Array{Float64,2}       = zeros(Float64, I, 2) # backward " "
    dV_upwind::Array{Float64,2} = zeros(Float64, I, 2) # value function difference matrix post upwind
    c::Array{Float64,2}         = similar(dVf) # flow consumption across state space
    u::Array{Float64,2}         = similar(dVf) # flow utility across state space
    g::Array{Float64,2}         = similar(dVf) # distribution vector
    r::Float64                  = r0  # initial interest rate guess
    KD::Float64 = (((α) / (r + δ)) ^ (1 / (1 - α))) * z_avg # capital demand
    w::Float64  = (1 - α) * (KD ^ α) * ((z_avg) ^ (-α))     # wage rate
    v0::Array{Float64,2} = similar(dVf)                     # initial value function guess
    v0[:, 1] = (w * μ * (1 - z[1]) .+ r .* a) .^ (1 - γ)/(1 - γ)/ρ # add guesses
    v0[:, 2] = (w * (1 - τ) * z[2] .+ r .* a) .^ (1 - γ)/(1 - γ)/ρ

    If::Array{Float64,2} = similar(dVf) # 1 if forward difference, 0 if not
    Ib::Array{Float64,2} = similar(dVf) # 1 if backward difference, 0 if not
    I0::Array{Float64,2} = similar(dVf)  # neither forward nor backward difference
    A::SparseMatrixCSC{Float64,Int64} = spzeros(Float64, 2*I, 2*I)

    # Iterate to find steady state interest rate
    for ir = 1:Ir

        # Get capital demand and wages
        KD = (((α) / (r + δ)) ^ (1 / (1 - α))) * z_avg
        w  = (1 - α) * (KD ^ α) * ((z_avg) ^ (-α))
        v  = v0

        # Solve for value function given interest rate
        for n = 1:maxit_HJB
            # Grab guess
            V = v

            # Compute forward difference
            dVf[1:I - 1, :] = (V[2:I, :] - V[1:(I - 1), :])/da
            dVf[I, :] = (w * ((1 - τ) * z .+ μ * (1 .- z)) .+ r .* amax).^(-γ)

            # Compute backward difference
            dVb[2:I, :] = (V[2:I, :] - V[1:(I - 1), :])/da
            dVb[1, :] = (w * ((1 - τ) * z .+ μ * (1 .- z)) .+ r .* amin).^(-γ)

            # Compute consumption and savings with forward difference
            cf  = dVf.^(-1/γ)
            ssf = w*((1 - τ) * zz + μ * (1 .- zz)) + r .* aa - cf

            # Compute consumption and savings with backward difference
            cb = dVb.^(-1/γ)
            ssb = w*((1 - τ) * zz + μ * (1 .- zz)) + r .* aa - cb

            # Compute consumption and derivative of value function for no drift
            c0::Array{Float64,2}  = w*((1 - τ) * zz .+ μ * (1 .- zz)) + r.*aa
            dV0::Array{Float64,2} = c0.^(-γ)

            # Compute upwind differences and matrix
            dV_upwind, If, Ib, I0 = upwind_value_function(dVf, dVb, dV0, ssf, ssb)
            c = dV_upwind .^ (-1 / γ)
            u = c .^ (1 - γ) / (1 - γ)
            savingsSS = w * ((1 - τ) * zz + μ * (1 .- zz)) + r .* aa .- c

            # Create A matrix and solve for hjb
            A = upwind_matrix(A_switch, ssf, ssb, da, I, J)
            V_stacked = solve_hjb(A, ρ, Δ, vec(u), vec(V))::Vector{Float64}
            V = reshape(V_stacked, I, J)

            # Update values and check convergence
            Vchange = V - v
            v = V

            if maximum(abs.(Vchange)) < crit_HJB
                break
            end
        end
        v0 = v

        # Solve for stationary distribution
        gg = solve_kfe(A, [da; dz], [I; J])::Vector{Float64}
        g = reshape(gg, I, J)

        # Compute capital supply, check market-clearing
        KS = Float64(g[:,1]' * a * da + g[:,2]' * a * da)
        S = Float64(KS - KD)

        # Update interest rate
        if S > crit_S
            # Excess Supply (for capital)
            rmax = r
            r = 0.5 * (r + rmin)
        elseif S < -crit_S
            # Excess Demand
            rmin = r
            r = 0.5 * (r + rmax)
        elseif abs(S) < crit_S
            # Set steady state values
            # Must maintain type of declared SteadyState parameter type
            m[:V_ss]  = vec(v)
            m[:gg_ss] = vec(g)[1:2*I - 1]
            m[:log_aggregate_tfp_ss] = 0.
            m[:K_ss]  = KS
            m[:r_ss]  = r
            m[:w_ss]  = w
            m[:Y_ss]  = (m[:K_ss]^α) * (z_avg ^ (1 - α))
            m[:C_ss]  = sum(c[:] .* g[:] * da)
            m[:I_ss]  = δ * m[:K_ss]
            m[:If_ss] = If
            m[:Ib_ss] = Ib
            m[:I0_ss] = I0
            break
        end
    end
    return m
end


function model_settings!(m::KrusellSmithCT)
    default_settings!(m)

    # State space grid
    m <= Setting(:I, 100, "Dimension of wealth grid")
    m <= Setting(:amin, 0., "Minimium wealth")
    m <= Setting(:amax, 100., "Maximum wealth")
    m <= Setting(:a, collect(range(0., stop=100., length=100)), "Wealth grid")
    m <= Setting(:da, 100. / (100 - 1), "Size of partitions in wealth grid")
    m <= Setting(:aa, [get_setting(m, :a) get_setting(m, :a)],
                 "Repetition of wealth grid across income grid dimension")
    m <= Setting(:aaa, reshape(get_setting(m, :aa), 2 * 100, 1), "Stacked repetition of wealth grid")
    m <= Setting(:J, 2, "Dimension of income grid")
    m <= Setting(:z, [0., 1.], "Income grid")
    m <= Setting(:dz, 1., "Size of partitions in income grid")
    m <= Setting(:zz, ones(100, 1) * get_setting(m, :z)', "Repetition of income grid across wealth dimension")
    m <= Setting(:zzz, reshape(get_setting(m, :zz), 2 * 100, 1), "Stacked repetition of income grid")
    endo = m.endogenous_states

    # Number of variables
    m <= Setting(:n_jump_vars, length(endo[:value_function]), "Number of jump variables")
    m <= Setting(:n_state_vars, length(endo[:distribution]) + n_shocks_exogenous(m),
                 "Number of state variables being reduced")
    m <= Setting(:n_state_vars_unreduce, 0, "Number of state variables not being reduced")

    # Approximation of steady state
    m <= Setting(:rmin, .0001, "Lower bound for steady state interest rate during bisection search")
    m <= Setting(:rmax, m[:ρ].value, "Upper bound for steady state interest rate during bisection search")
    m <= Setting(:r0, .005, "Initial guess for steady state interest rate")
    m <= Setting(:maxit_HJB, 100, "Maximum iterations on steady state HJB")
    m <= Setting(:crit_HJB, 1e-6, "Error criterion for convergence of steady state HJB")
    m <= Setting(:Δ_HJB, 1e4, "Update size for implicit scheme for steady state HJB")
    m <= Setting(:Ir, 100, "Maximum iterations for steady state interest rate")
    m <= Setting(:crit_S, 1e-5, "Error criterion for convergence of steady state interest rate")

    # Reduction
    m <= Setting(:reduce_state_vars, true, "Reduce state variables or not")
    m <= Setting(:reduce_v, true, "Reduce value function or not")
    m <= Setting(:krylov_dim, 5, "Krylov reduction dimension")
    m <= Setting(:F, identity, "Function applied during block arnoldi deflation")
    m <= Setting(:n_knots, 12, "Number of knot points")
    m <= Setting(:c_power, 7, "Amount of bending of knot point locations to make them nonuniform")

    # Create knot points for spline reductoin
    knots = collect(range(get_setting(m, :amin),
                          stop=get_setting(m, :amax), length=get_setting(m, :n_knots) .- 1))
    knots = (get_setting(m, :amax) - get_setting(m, :amin))/(2^get_setting(m, :c_power) - 1) *
        ((knots .- get_setting(m, :amin)) / (get_setting(m, :amax) .- get_setting(m, :amin)) .+ 1) .^
        get_setting(m, :c_power) .+ get_setting(m, :amin) .- (get_setting(m, :amax) .-
                                   get_setting(m, :amin))/(2^get_setting(m, :c_power) - 1)

    # Reduction continued
    m <= Setting(:knots_dict, Dict(1 => collect(range(get_setting(m, :amin),
                                   stop=get_setting(m, :amax), length=(12 - 1)))),
                 "Location of knot points for each dimension for value function reduction")
    m <= Setting(:spline_grid, get_setting(m, :a), "Grids for spline basis")
    m <= Setting(:n_prior, 1,
         "Number of dimensions approximated by spline basis that were not used to compute the basis matrix")
    m <= Setting(:n_post, 2, "Number of dimensions that need to be approximated by spline basis")

    # Transition time step
    m <= Setting(:dt, 1. / 2, "Size of time step for transition equation")

    # Track last lag or not?
    m <= Setting(:track_lag, false, "Add first lag when constructing measurement equation")

    # Sampling Method
    m <= Setting(:sampling_method, :SMC, "Set sampling method to SMC")

    # SMC settings
    m <= Setting(:n_particles, 3_000)
    m <= Setting(:n_Φ, 200)
    m <= Setting(:λ, 2.0)
    m <= Setting(:n_smc_blocks, 1)
    m <= Setting(:use_parallel_workers, true)
    # m <= Setting(:use_parallel_workers, false)
    m <= Setting(:step_size_smc, 0.5)
    m <= Setting(:n_mh_steps_smc, 1)
    m <= Setting(:resampler_smc, :polyalgo)
    m <= Setting(:target_accept, 0.25)
    m <= Setting(:mixture_proportion, .9)
    m <= Setting(:tempering_target, 0.95)
    m <= Setting(:resampling_threshold, .5)
    m <= Setting(:use_fixed_schedule, true)
    m <= Setting(:smc_iteration, 0)

    # Forecast
    m <= Setting(:use_population_forecast, true,
                 "Whether to use population forecasts as data")
    m <= Setting(:forecast_zlb_value, 0.13,
        "Value of the zero lower bound in forecast periods, if we choose to enforce it")

    #Simulating states
    m <= Setting(:state_simulation_freq, 2,
                 "How many states you want to simulate between states + 1")

end

"""
Overloading functions so as to properly grab information from CT Hank Models.
"""
n_states(m::KrusellSmithCT) = sum(map(i -> length(collect(m.endogenous_states)[i][2]), 1:length(keys(m.endogenous_states))))
n_shocks_expectational(m::KrusellSmithCT) = sum(map(i -> length(collect(m.expected_shocks)[i][2]), 1:length(keys(m.expected_shocks))))

#=
n_states(m::OneAssetHANK) = sum(map(i -> length(collect(m.endogenous_states)[i][2]), 1:length(keys(m.endogenous_states))))
n_shocks_expectational(m::OneAssetHANK) = sum(map(i -> length(collect(m.expected_shocks)[i][2]), 1:length(keys(m.expected_shocks))))
=#
