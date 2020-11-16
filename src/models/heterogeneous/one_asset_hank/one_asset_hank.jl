"""
```
OneAssetHANK{T} <: AbstractCTModel{T}
```

The `OneAssetHANK` type defines the structure of the one-asset HANK model
originally written in MATLAB by Ben Moll and SeHyoun Ahn at
https://sehyoun.com/EXAMPLE_one_asset_HANK_web.html

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

* `spec::String`: The model specification identifier, \"one_asset_hank\", cached
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
mutable struct OneAssetHANK{T} <: AbstractCTModel{T}
    parameters::ParameterVector{T}                          # vector of all time-invariant model parameters
    steady_state::ParameterVector{T}                        # model steady-state values
    keys::OrderedDict{Symbol,Int}                           # human-readable names for all the model
                                                            # parameters and steady-states

    endogenous_states::OrderedDict{Symbol,Vector{Int}}      # these fields used to create matrices in the
    exogenous_shocks::OrderedDict{Symbol,Int}               # measurement & equilibrium condition equations
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

description(m::OneAssetHANK) = "Julia implementation of the one-asset HANK model
originally written in MATLAB by Ben Moll and SeHyoun Ahn at
https://sehyoun.com/EXAMPLE_one_asset_HANK_web.html: OneAssetHANK, $(m.subspec)"


"""
`init_model_indices!(m::OneAssetHANK)`

Arguments:
`m:: OneAssetHANK`: a model object

Description:
Initializes indices for all of `m`'s states, shocks, and equilibrium conditions.
"""
function init_model_indices!(m::OneAssetHANK)
    # Endogenous states
    endogenous_states = collect([:value_function, :inflation, :distribution,
                                 :monetary_policy, :w, :N, :C, :output, :B])

    # Exogenous shocks
    exogenous_shocks = collect([:mp_shock])

    # Expectations shocks
    expected_shocks = collect([:E_V, :E_π])

    # Equilibrium conditions
    equilibrium_conditions = collect([:value_function, :inflation, :distribution,
                                      :monetary_policy, :w, :N, :C, :Y, :B])

    # Additional states added after solving model
    # Lagged states and observables measurement error
    endogenous_states_augmented = []

    # Observables
    observables = keys(m.observable_mappings)

    # Pseudo-observables
    pseudo_observables = keys(m.pseudo_observable_mappings)

    # Assign dimensions
    m.endogenous_states[:value_function]  = collect(1:200)      # size of state space grid
    m.endogenous_states[:inflation]       = [201]               # inflation
    m.endogenous_states[:distribution]    = collect(202:400)    # size of grid minus 1 for ...
    m.endogenous_states[:monetary_policy] = [401]               # monetary policy shocks
    m.endogenous_states[:w]               = [402]               # static conditions from market-clearing
    m.endogenous_states[:N]               = [403]
    m.endogenous_states[:C]               = [404]
    m.endogenous_states[:output]          = [405]
    m.endogenous_states[:B]               = [406]
    m.expected_shocks[:E_V]               = m.endogenous_states[:value_function]
    m.expected_shocks[:E_π]               = m.endogenous_states[:inflation]

    for (i,k) in enumerate(endogenous_states_augmented)
        m.endogenous_states_augmented[k] = i + length(endogenous_states)
    end
    for (i,k) in enumerate(observables);          m.observables[k]          = i end
    for (i,k) in enumerate(exogenous_shocks);     m.exogenous_shocks[k]     = i end
    for (i,k) in enumerate(pseudo_observables);   m.pseudo_observables[k]   = i end
end


function OneAssetHANK(subspec::String="ss0";
                       custom_settings::Dict{Symbol, Setting} = Dict{Symbol, Setting}(),
                       testing = false)

    # Model-specific specifications
    spec               = split(basename(@__FILE__),'.')[1]
    subspec            = subspec
    settings           = Dict{Symbol,Setting}()
    test_settings      = Dict{Symbol,Setting}()
    rng                = MersenneTwister(0)

    # Initialize empty model
    m = OneAssetHANK{Float64}(
            # Model parameters and steady state values # before 2nd arg was Vector{Float64}()
            Vector{AbstractParameter{Float64}}(), Vector{Float64}(),
            OrderedDict{Symbol,Int}(),

            # Model indices
            OrderedDict{Symbol,Vector{Int}}(), OrderedDict{Symbol,Int}(),
            OrderedDict{Symbol,Vector{Int}}(), OrderedDict{Symbol,Vector{Int}}(),
            OrderedDict{Symbol,Int}(),         OrderedDict{Symbol,Int}(),
            OrderedDict{Symbol,Int}(),

            spec,
            subspec,
            settings,
            test_settings,
            rng,
            testing,
            OrderedDict{Symbol,Observable}(),
            OrderedDict{Symbol,PseudoObservable}())

    # Initialize parameters
    init_parameters!(m)
    init_model_indices!(m)

    # Set settings
    model_settings!(m)
    default_test_settings!(m)

    for custom_setting in values(custom_settings)
        m <= custom_setting
    end

    # Test to see if necessary values are defined, throw error otherwise
    try
        @assert typeof(get_setting(m, :reduce_state_vars)) == Bool
        @assert typeof(get_setting(m, :reduce_v)) == Bool
    catch
        error("Need to specify in settings whether to perform distribution
               and/or value function reduction")
    end

    if get_setting(m, :reduce_state_vars) || get_setting(m, :reduce_v)
        try
            @assert get_setting(m, :n_jump_vars) > 0
            @assert get_setting(m, :n_state_vars) > 0
            @assert get_setting(m, :n_state_vars_unreduce) >= 0
        catch
            error("Need to enter values in settings for n_jump_vars, n_state_vars,
                   and/or n_state_vars_unreduce in Settings. These indicate the
                   number of jump and aggregate state variables.")
        end
    end

    if get_setting(m, :reduce_state_vars)
        try
            @assert get_setting(m, :krylov_dim) > 0
        catch
            error("Need to specify in settings desired dimension of
                   Krylov subspace, krylov_dim.")
        end
    end

    if get_setting(m, :reduce_v)
        try
            @assert get_setting(m, :n_prior) > 0 && typeof(get_setting(m, :n_prior)) == Int
            @assert get_setting(m, :n_post) > 0 && typeof(get_setting(m, :n_post)) == Int
            @assert in(:knots_dict, keys(m.settings))
            @assert in(:spline_grid, keys(m.settings))
        catch
            error("Need to specify in settings n_prior, n_post, knots_dict
                   and/or spline_grid for value function reduction.
                   See reduction.jl for their uses.")
        end
    end

    # Compute steady state
    steadystate!(m)
    return m
end

"""
```
init_parameters!(m::OneAssetHANK)
```

Initializes the model's parameters, as well as empty values for the steady-state
parameters (in preparation for `steadystate!(m)` being called to initialize
those).
"""
function init_parameters!(m::OneAssetHANK)
    # Initialize parameters
    m <= parameter(:coefrra, 1.0, fixed=true,
                   description="Relative risk aversion",
                   tex_label="coefrra")

    m <= parameter(:frisch, 0.5, fixed=true,
                   description="Frisch elasticity",
                   tex_label="frisch")

    m <= parameter(:meanlabeff, 3.0, fixed=true,
                   description="meanlabeff: so that at h=1/3 output will be approximately = 1",
                   tex_label="meanlabeff")

    m <= parameter(:maxhours, 1.0, fixed=true,
                   description="maxhours:...",
                   tex_label="maxhours")

    # Production
    m <= parameter(:ceselast, 10., fixed=true,
                   description="CES elasticity",
                   tex_label="ceselast")

    m <= parameter(:priceadjust, 100., fixed=true,
                   description="priceadjust...",
                   tex_label="priceadjust")

    # Policy parameters
    m <= parameter(:taylor_inflation, 1.25, fixed=true,
                   description="Taylor rule coefficient on inflation",
                   tex_label="taylor_inflation")

    m <= parameter(:taylor_outputgap, 0., fixed=true,
                   description="Taylor rule coefficient on output",
                   tex_label="taylor_outputgap")

    m <= parameter(:labtax, 0.2, fixed=true,
                   description="Marginal tax rate on labor income",
                   tex_label="labtax")

    m <= parameter(:govbondtarget, 6., fixed=true,
                   description="govbondtarget: multiple of quarterly GDP",
                   tex_label="govbondtarget")

    m <= parameter(:lumptransferpc, 0.06, fixed=true,
                   description="lumptransferepc: 6% of quarterly GDP in steady state",
                   tex_label="lumptransferepc")

    m <= parameter(:govbcrule_fixnomB, 0., fixed=true,
                   description="govbcrule_fixnomB...",
                   tex_label="govbvrule_fixnomB")

    m <= parameter(:labdisutil, m[:meanlabeff] / ((0.75 ^(-m[:coefrra])) *
                                                  ((1. / 3.)^(1/m[:frisch]))), fixed=true,
                   description="Coefficient of labor disutility",
                   tex_label="labdisutil")

    # Aggregate shocks
    m <= parameter(:σ_MP, sqrt(0.05), (1e-20, 1e5), (1e-20, 1e5), ModelConstructors.Exponential(),
                   RootInverseGamma(4, .4), fixed=false,
                   description="Volatility of monetary policy shocks", tex_label="σ_MP")

    m <= parameter(:θ_MP, 0.25, fixed=true,
                   description="Rate of mean reversion in monetary policy shocks",
                   tex_label="θ_MP")

    m <= SteadyStateParameterGrid(:V_ss, Vector{Float64}(),
                                  description = "Stacked steady-state value function")

    m <= SteadyStateParameter(:inflation_ss, NaN, description = "Steady-state rate of inflation")
    m <= SteadyStateParameterGrid(:g_ss, Vector{Float64}(),
                                  description = "Stacked steady-state distribution")
    m <= SteadyStateParameter(:r_ss, NaN, description = "Steady-state real interest rate")
    m <= SteadyStateParameterGrid(:u_ss, Vector{Float64}(), description = "Steady-state flow utility")
    m <= SteadyStateParameterGrid(:c_ss, Vector{Float64}(),
                                  description = "Steady-state individual consumption rate")
    m <= SteadyStateParameterGrid(:h_ss, Vector{Float64}(),
                                  description = "Steady-state individual labor hours")
    m <= SteadyStateParameterGrid(:s_ss, Vector{Float64}(),
                                  description = "Steady-state individual savings rate")
    m <= SteadyStateParameter(:rnom_ss, NaN, description = "Steady-state nominal interest rate")
    m <= SteadyStateParameter(:B_ss, NaN, description = "Steady-state bond supply")
    m <= SteadyStateParameter(:N_ss, NaN, description = "Steady-state labor supply")
    m <= SteadyStateParameter(:Y_ss, NaN, description = "Steady-state output")
    m <= SteadyStateParameter(:labor_share_ss, NaN, description = "Frictionless labor share of income")
    m <= SteadyStateParameter(:w_ss, NaN, description = "Steady-state wage rate")
    m <= SteadyStateParameter(:profit_ss, NaN, description = "Steady-state profits")
    m <= SteadyStateParameter(:C_ss, NaN, description = "Steady-state aggregate consumption")
    m <= SteadyStateParameter(:T_ss, NaN, description = "Steady-state total taxes")
    m <= SteadyStateParameter(:G_ss, NaN, description = "Steady-state government spending")
    m <= SteadyStateParameter(:ρ_ss, NaN, description = "Steady-state discount rate")
end

"""
```
steadystate!(m::OneAssetHANK)
```

Calculates the model's steady-state values. `steadystate!(m)` must be called whenever
the parameters of `m` are updated.
"""
function steadystate!(m::OneAssetHANK)

    # Read in parameters
    coefrra        = m[:coefrra].value
    frisch         = m[:frisch].value
    meanlabeff     = m[:meanlabeff].value
    maxhours       = m[:maxhours].value
    ceselast       = m[:ceselast].value
    labtax         = m[:labtax].value
    govbondtarget  = m[:govbondtarget].value
    labdisutil     = m[:labdisutil].value
    lumptransferpc = m[:lumptransferpc].value

    # Read in grids
    I                = get_setting(m, :I)
    J                = get_setting(m, :J)
    a                = get_setting(m, :a)
    g_z              = get_setting(m, :g_z)
    zz               = get_setting(m, :zz)
    ymarkov_combined = get_setting(m, :ymarkov_combined)

    # Set necessary variables
    aa   = repeat(a, 1, J)
    amax = maximum(a)
    amin = minimum(a)

    # Read in initial rates
    iterate_r = get_setting(m, :iterate_r)
    r         = get_setting(m, :r0)
    rmin      = get_setting(m, :rmin)
    rmax      = get_setting(m, :rmax)
    iterate_ρ = get_setting(m, :iterate_ρ)
    ρ         = get_setting(m, :ρ0)
    ρmin      = get_setting(m, :ρmin)
    ρmax      = get_setting(m, :ρmax)

    # Read in approximation parameters
    Ir          = get_setting(m, :Ir)
    maxit_HJB   = get_setting(m, :maxit_HJB)
    tol_HJB     = get_setting(m, :tol_HJB)
    Δ_HJB       = get_setting(m, :Δ_HJB)
    maxit_kfe   = get_setting(m, :maxit_kfe)
    tol_kfe     = get_setting(m, :tol_kfe)
    Δ_kfe       = get_setting(m, :Δ_kfe)
    niter_hours = get_setting(m, :niter_hours)
    crit_S      = get_setting(m, :crit_S)

    # Initializing equilibrium objects
    labor_share_ss = (ceselast - 1) / ceselast
    w       = w_ss = labor_share_ss

    # compute initial guesses at steady state values given zz, labor_share_ss, etc.
    N_ss, Y_ss, B_ss, profit_ss, profshare, lumptransfer =
        calculate_ss_equil_vars(zz, labor_share_ss, meanlabeff, lumptransferpc, govbondtarget)

    # Initialize matrices for finite differences
    Vaf = Array{ComplexF64}(undef, I, J)
    Vab = Array{ComplexF64}(undef, I, J)

    cf  = Array{ComplexF64}(undef, I, J) # forward consumption difference
    hf  = Array{ComplexF64}(undef, I, J) # forward hours difference
    sf  = Array{ComplexF64}(undef, I, J) # forward saving difference
    cb  = Array{ComplexF64}(undef, I, J) # backward consumption difference
    hb  = Array{ComplexF64}(undef, I, J) # backward hours difference
    sb  = Array{ComplexF64}(undef, I, J) # backward saving difference
    c0  = Array{ComplexF64}(undef, I, J)
    A   = Array{ComplexF64}(undef, I*J, I*J)

    # Aswitch*v = \lambda_z(v(a,z') - v(a,z))
    # Captures expected change in value function due to jumps in z
    Aswitch = kron(ymarkov_combined, SparseMatrixCSC{ComplexF64}(LinearAlgebra.I, I, I))

    # Initialize steady state variables
    V  = Array{ComplexF64}(undef, I, J) # value function
    u  = Array{ComplexF64}(undef, I, J) # flow utility across state space
    s  = Array{ComplexF64}(undef, I, J) # savings across state space
    c  = Array{ComplexF64}(undef, I, J) # flow consumption
    h  = Array{ComplexF64}(undef, I, J) # flow hours of labor
    h0 = Array{ComplexF64}(undef, I, J) # guess of what h will be

    # Creates functions for computing flow utility, income earned, and labor done given
    # CRRA + frisch elasticity style labor disutility
    util, income, labor = construct_household_problem_functions(V, w, coefrra, frisch, labtax, labdisutil)

    # Setting up forward/backward difference grids for a. See helpers.jl
    daf, dab, azdelta = initialize_diff_grids(a, I, J)

    for ir = 1:Ir
        fill!(c,  complex(0.))
        fill!(h,  complex(1/3))
        fill!(h0, complex(1.))

        # Initial guess
        v = similar(h)

        for i in eachindex(h)
            inc  = income(h[i], zz[i], profshare[i], lumptransfer, r, aa[i]) # Get income
            v[i] = util(inc, h[i]) / ρ                                       # Value function guess
        end

        # Iterate HJB
        for ihjb = 1:maxit_HJB
            V = v
            Vaf, Vab, cf, hf, cb, hb = construct_initial_diff_matrices(V, Vaf, Vab,
                                                                       income, labor, h, h0, zz,
                                                                       profshare, lumptransfer,
                                                                       amax, amin,
                                                                       coefrra, r, daf,
                                                                       dab, maxhours)

            # Iterative method to find consistent forward/backward/neutral difference matrices for c and h
            cf, hf, cb, hb, c0, h0 = hours_iteration(income, labor, zz, profshare,
                                                     lumptransfer, aa,
                                                     coefrra, r, cf, hf, cb, hb, c0, h0,
                                                     maxhours, niter_hours)

            for i in eachindex(h0)
                c0[i] = income(h0[i], zz[i], profshare[i], lumptransfer, r, aa[i])
                sf[i] = income(hf[i], zz[i], profshare[i], lumptransfer, r, aa[i]) - cf[i]
                sb[i] = income(hb[i], zz[i], profshare[i], lumptransfer, r, aa[i]) - cb[i]
            end

            for j=1:J
                Vaf[I, j] = cf[I, j] ^ (-coefrra) # Forward difference for value function w.r.t. wealth
                Vab[1, j] = cb[1, j] ^ (-coefrra) # Backward difference for value function w.r.t. wealth
            end

            # See helpers.jl: this is a custom upwind function, not one from solve_hjb.jl
            V, A, u, h, c, s  = upwind(ρ, V, util, Aswitch, cf, cb, c0, hf, hb, h0, sf, sb, Vaf, Vab,
                                    daf, dab; Δ_HJB = Δ_HJB)

            # Check for convergence
            Vchange = V - v
            v = V

            err_HJB = maximum(abs.(Vchange))
            if err_HJB < tol_HJB
                break
            end
        end

        # Create initial guess for g0
        g0 = zeros(ComplexF64, I, J)
        # Assign stationary income distribution weight at a = 0, zero elsewhere
        g0[iszero.(a), :] = g_z
        # g_z is marginal distribution, so re-weight by some multiplier of Lebesgue measure
        g0 = g0 ./ reshape(azdelta, I, J)
        # Solve for distribution
        g = solve_kfe(A, g0, spdiagm(0 => azdelta),
                       maxit_kfe = maxit_kfe, tol_kfe = tol_kfe,
                       Δ_kfe = Δ_kfe)

        # Back out static conditions/steady state values given our value function and distribution
        N_ss, Y_ss, B_ss, profit_ss, profshare, lumptransfer, bond_err =
            calculate_ss_equil_vars(zz, h, g, azdelta, aa, labor_share_ss, meanlabeff,
                                lumptransferpc, govbondtarget)

        # Check bond market for market clearing
        r, rmin, rmax, ρ, ρmin, ρmax, clear_cond =
            check_bond_market_clearing(bond_err, crit_S, r, rmin, rmax, ρ, ρmin, ρmax,
                                       iterate_r, iterate_ρ)

        if clear_cond
            # Set steady state values
            m[:V_ss]           = real(vec(V))
            m[:inflation_ss]   = 0.
            m[:g_ss]           = real(g)
            m[:r_ss]           = r
            m[:u_ss]           = real(vec(u))
            m[:c_ss]           = real(vec(c))
            m[:h_ss]           = real(vec(h))
            m[:s_ss]           = real(vec(s))
            m[:rnom_ss]        = m[:r_ss].value + m[:inflation_ss].value
            m[:B_ss]           = sum(m[:g_ss].value .* vec(aa) .* azdelta)
            m[:N_ss]           = real(sum(vec(zz) .* m[:h_ss].value .* m[:g_ss].value .* azdelta))
            m[:Y_ss]           =  m[:N_ss]
            m[:labor_share_ss] = (ceselast - 1) / ceselast
            m[:w_ss]           =  m[:labor_share_ss].value
            m[:profit_ss]      = (1 - m[:labor_share_ss].value) .* m[:Y_ss].value
            m[:C_ss]           = sum(m[:c_ss].value .* m[:g_ss].value .* azdelta)
            m[:T_ss]           = real(lumptransfer)
            m[:ρ_ss]           = ρ
            m[:G_ss]           = labtax * m[:w_ss].value * m[:N_ss].value -
                                 m[:T_ss].value - m[:r_ss].value * m[:B_ss].value
            break
        end
    end
    return m
end

function model_settings!(m::OneAssetHANK)
    default_settings!(m)

    # Steady state r iteration
    m <= Setting(:r0, 0.005, "Initial guess for real interest rate")
    m <= Setting(:rmax, 0.08, "Maximum guess for real interest rate")
    m <= Setting(:rmin, 0.001, "Minimum guess for real interest rate")

    # Steady state ρ iteration
    m <= Setting(:ρ0, 0.02, "Initial guess for discount rate")
    m <= Setting(:ρmax, 0.05, "Maximum guess for discount rate")
    m <= Setting(:ρmin, 0.005, "Minimum guess for discount rate")

    # State space grid
    m <= Setting(:I, 100, "Number of grid points")
    m <= Setting(:J, 2, "Number of income states")
    m <= Setting(:amin, 0., "Minimum asset grid value")
    m <= Setting(:amax, 40., "Maximum asset grid value")
    m <= Setting(:agridparam, 1, "Bending coefficient: 1 for linear")
    m <= Setting(:a, construct_asset_grid(get_setting(m,:I), get_setting(m, :agridparam),
                        get_setting(m, :amin), get_setting(m,:amax)), "Asset grid")
    m <= Setting(:ygrid_combined, [0.2, 1])
    m <= Setting(:ymarkov_combined, [-0.5 0.5; 0.0376 -0.0376], "Markov transition parameters")
    m <= Setting(:g_z, compute_stationary_income_distribution(get_setting(m, :ymarkov_combined),
                       get_setting(m, :J)), "Stationary income distribution")
    m <= Setting(:zz, construct_labor_income_grid(get_setting(m, :ygrid_combined),
                       get_setting(m, :g_z), m[:meanlabeff].value, get_setting(m, :I)),
                       "Labor income grid repeated across asset dimension")
    m <= Setting(:z, get_setting(m, :zz)[1, :], "Labor income grid")

    # Number of variables
    m <= Setting(:n_jump_vars, get_setting(m, :I) * get_setting(m, :J) + 1, "Number of jump variables")
    m <= Setting(:n_state_vars, get_setting(m, :I) * get_setting(m, :J), "Number of state variables")
    m <= Setting(:n_state_vars_unreduce, 0, "Number of state variables not being reduced")
    # R: Inserted from old code; not certain necessary
    m <= Setting(:n_static_relations, 5, "Number of static relations: bond-market clearing, labor market
                                          clearing, consumption, output, total assets")
    m <= Setting(:n_vars, get_setting(m, :n_jump_vars) + get_setting(m, :n_state_vars)
                 + get_setting(m, :n_static_relations),
                 "Number of variables, total")

    # Steady state approximation
    m <= Setting(:Ir, 100, "Max number of iterations for steady state")
    m <= Setting(:maxit_HJB, 500, "Max number of iterations for HJB")
    m <= Setting(:tol_HJB, 1e-8, "Tolerance for HJB error")
    m <= Setting(:Δ_HJB, 1e6, "Multiplier for implicit upwind scheme of HJB")
    m <= Setting(:maxit_kfe, 1000, "Max number of iterations for kfe")
    m <= Setting(:tol_kfe, 1e-12, "Tolerance for kfe error")
    m <= Setting(:Δ_kfe, 1e6, "Multiplier for implicit upwind scheme of kfe")
    m <= Setting(:niter_hours, 10, "Max number of iterations for finding hours worked")
    m <= Setting(:iterate_r, false, "Iterate on real interest rate or not")
    m <= Setting(:iterate_ρ, true, "Iterate on discount rate or not")
    m <= Setting(:crit_S, 1e-5, "Tolerance for error in bond markets")

    # Reduction
    m <= Setting(:n_knots, 12, "Number of knot points")
    m <= Setting(:c_power, 1, "Amount of bending of knot point locations to make them nonuniform")
    m <= Setting(:n_post, length(get_setting(m, :zz)[1,:]),
                 "Number of dimensions that need to be approximated by spline basis")
    m <= Setting(:n_prior, 1,
                 "Number of dimensions approximated by spline basis that
                  were not used to compute the basis matrix")

    knots = collect(range(get_setting(m, :amin), stop=get_setting(m, :amax), length=get_setting(m, :n_knots)-1))
    knots = (get_setting(m, :amax) - get_setting(m, :amin)) /
            (2^get_setting(m, :c_power)-1) * ((knots .- get_setting(m, :amin))
                                              / (get_setting(m, :amax) -
                                                 get_setting(m, :amin)) .+ 1) .^ get_setting(m, :c_power) .+
                                                 get_setting(m, :amin) .- (get_setting(m, :amax) -
                                                                          get_setting(m, :amin)) /
                                                                          (2 ^ get_setting(m, :c_power) - 1)
    m <= Setting(:knots_dict, Dict(1 => knots),
                 "Location of knot points for each dimension for value function reduction")
    m <= Setting(:krylov_dim, 20, "Krylov reduction dimension")
    m <= Setting(:reduce_state_vars, true, "Reduce state variables or not") #turned off reduction
    m <= Setting(:reduce_v, true, "Reduce value function or not") #turned off reduction
    m <= Setting(:spline_grid, get_setting(m, :a), "Grid of knot points for spline basis")

    # Sampling method
    m <= Setting(:sampling_method, :SMC, "Set sampling method to SMC")

    # SMC Settings
    m <= Setting(:n_particles, 3000)
    m <= Setting(:n_Φ, 200)
    m <= Setting(:λ, 2.0)
    m <= Setting(:n_smc_blocks, 1)
    m <= Setting(:use_parallel_workers, true)
    m <= Setting(:step_size_smc, 0.5)
    m <= Setting(:n_mh_steps_smc, 1)
    m <= Setting(:resampler_smc, :polyalgo)
    m <= Setting(:target_accept, 0.25)
    m <= Setting(:mixture_proportion, .9)
    m <= Setting(:tempering_target, 0.95)
    m <= Setting(:resampling_threshold, .5)
    m <= Setting(:use_fixed_schedule, true)
    m <= Setting(:smc_iteration, 0)

    # Track lag
    m <= Setting(:track_lag, false, "Add first lag when constructing measurement equation")

    # Forecast
    m <= Setting(:use_population_forecast, true,
                 "Whether to use population forecasts as data")
    m <= Setting(:forecast_zlb_value, 0.13,
        "Value of the zero lower bound in forecast periods, if we choose to enforce it")

    # Simulating states
    m <= Setting(:state_simulation_freq, 2,
                 "How many states you want to simulate between states + 1")
end
