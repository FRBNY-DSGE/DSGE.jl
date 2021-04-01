"""
```
Model1002{T} <: AbstractRepModel{T}
```

The `Model1002` type defines the structure of Model1002 (same as
Model997 but uses longrate without adjusting for term premia.)

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

* `spec::String`: The model specification identifier, \"m1002\", cached here for
  filepath computation.

* `subspec::String`: The model subspecification number, indicating that
  some parameters from the original model spec (\"ss10\") are initialized
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
  stores data sources, series mnemonics, and transformations to/from model
  units. DSGE.jl will fetch data from the Federal Reserve Bank of St. Louis's
  FRED database; all other data must be downloaded by the user. See `load_data`
  and `Observable` for further details.

* `pseudo_observable_mappings::OrderedDict{Symbol,PseudoObservable}`: A
  dictionary that stores names and transformations to/from model units. See
  `PseudoObservable` for further details.
"""
mutable struct Model1002{T} <: AbstractRepModel{T}
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

description(m::Model1002) = "New York Fed DSGE Model m1002, $(m.subspec)."

"""
`init_model_indices!(m::Model1002)`

Arguments:
`m:: Model1002`: a model object

Description:
Initializes indices for all of `m`'s states, shocks, and equilibrium conditions.
"""
function init_model_indices!(m::Model1002)
    # Endogenous states
    endogenous_states = [[
        :y_t, :c_t, :i_t, :qk_t, :k_t, :kbar_t, :u_t, :rk_t, :Rktil_t, :n_t, :mc_t,
        :π_t, :μ_ω_t, :w_t, :L_t, :R_t, :g_t, :b_t, :μ_t, :z_t, :λ_f_t, :λ_f_t1,
        :λ_w_t, :λ_w_t1, :rm_t, :σ_ω_t, :μ_e_t, :γ_t, :π_star_t,
        :Ec_t, :Eqk_t, :Ei_t, :Eπ_t, :EL_t, :Erk_t, :Ew_t, :ERktil_t, :ERktil_f_t,
        :y_f_t, :c_f_t, :i_f_t, :qk_f_t, :k_f_t,
        :kbar_f_t, :u_f_t, :rk_f_t, :w_f_t, :L_f_t, :r_f_t, :Ec_f_t, :Eqk_f_t, :Ei_f_t, :EL_f_t,
        :ztil_t, :π_t1, :π_t2, :π_a_t, :R_t1, :zp_t, :Ez_t,
        :Rktil_f_t, :n_f_t];
        [Symbol("rm_tl$i") for i = 1:n_mon_anticipated_shocks(m)]]
    for (key, val) in get_setting(m, :antshocks)
        endogenous_states = vcat(endogenous_states, [Symbol(key, "_tl$i") for i = 1:val])
    end

    # Exogenous shocks
    exogenous_shocks = [[
        :g_sh, :b_sh, :μ_sh, :ztil_sh, :λ_f_sh, :λ_w_sh, :rm_sh, :σ_ω_sh, :μ_e_sh,
        :γ_sh, :π_star_sh, :zp_sh, :lr_sh, :tfp_sh, :gdpdef_sh, :corepce_sh, :gdp_sh, :gdi_sh];
        [Symbol("rm_shl$i") for i = 1:n_mon_anticipated_shocks(m)]]
    for (key, val) in get_setting(m, :antshocks)
        exogenous_shocks = vcat(exogenous_shocks, [Symbol(key, "_shl$i") for i = 1:val])
    end

    # Expectations shocks
    expected_shocks = [
        :Ec_sh, :Eqk_sh, :Ei_sh, :Eπ_sh, :EL_sh, :Erk_sh, :Ew_sh, :ERktil_sh, :Ec_f_sh,
        :Eqk_f_sh, :Ei_f_sh, :EL_f_sh, :ERktil_f_sh]

    # Equilibrium conditions
    equilibrium_conditions = [[
        :eq_euler, :eq_inv, :eq_capval, :eq_spread, :eq_nevol, :eq_output, :eq_caputl, :eq_capsrv, :eq_capev,
        :eq_mkupp, :eq_phlps, :eq_caprnt, :eq_msub, :eq_wage, :eq_mp, :eq_res, :eq_g, :eq_b, :eq_μ, :eq_z,
        :eq_λ_f, :eq_λ_w, :eq_rm, :eq_σ_ω, :eq_μ_e, :eq_γ, :eq_λ_f1, :eq_λ_w1, :eq_Ec,
        :eq_Eqk, :eq_Ei, :eq_Eπ, :eq_EL, :eq_Erk, :eq_Ew, :eq_ERktil, :eq_euler_f, :eq_inv_f,
        :eq_capval_f, :eq_output_f, :eq_caputl_f, :eq_capsrv_f, :eq_capev_f, :eq_mkupp_f,
        :eq_caprnt_f, :eq_msub_f, :eq_res_f, :eq_Ec_f, :eq_Eqk_f, :eq_Ei_f, :eq_EL_f,
        :eq_ztil, :eq_π_star, :eq_π1, :eq_π2, :eq_π_a, :eq_Rt1, :eq_zp, :eq_Ez, :eq_spread_f,:eq_nevol_f,  :eq_ERktil_f];
        [Symbol("eq_rml$i") for i=1:n_mon_anticipated_shocks(m)]]
    for (key, val) in get_setting(m, :antshocks)
        equilibrium_conditions = vcat(equilibrium_conditions, [Symbol("eq_", key, "l$i") for i = 1:val])
    end

    # Additional states added after solving model
    # Lagged states and observables measurement error
    endogenous_states_augmented = [
        :y_t1, :c_t1, :i_t1, :w_t1, :π_t1_dup, :L_t1, :u_t1, :Et_π_t, :e_lr_t, :e_tfp_t, :e_gdpdef_t,
        :e_corepce_t, :e_gdp_t, :e_gdi_t, :e_gdp_t1, :e_gdi_t1]
    if subspec(m) in ["ss13", "ss14", "ss15", "ss16", "ss17", "ss18", "ss19", "ss20"]
        push!(endogenous_states_augmented, :Sinf_t, :πtil_t, :πtil_t1)
    end
    if subspec(m) in ["ss14", "ss15", "ss16", "ss18", "ss19"]
        push!(endogenous_states_augmented, :e_tfp_t1)
    end

    # COVID-19 states, shocks, and equations
    if subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63"]
        push!(endogenous_states, :ziid_t)
        push!(equilibrium_conditions, :eq_ziid)
        push!(exogenous_shocks, :ziid_sh)
        push!(endogenous_states, :biidc_t)
        push!(equilibrium_conditions, :eq_biidc)
        push!(exogenous_shocks, :biidc_sh)
        push!(endogenous_states, :φ_t)
        push!(endogenous_states, :Eφ_t)
        push!(equilibrium_conditions, :eq_φ)
        push!(equilibrium_conditions, :eq_Eφ)
        push!(exogenous_shocks, :φ_sh)
    end

    # Add various pseudo-observables
    if get_setting(m, :add_laborproductivity_measurement)
        push!(endogenous_states_augmented, :cum_z_t)
        m <= Setting(:integrated_series, [:cum_z_t])
    end
    if get_setting(m, :add_nominalgdp_level)
        integ_series = [:cum_z_t, :cum_y_t, :cum_e_gdp_t, :cum_π_t]
        push!(endogenous_states_augmented, setdiff(integ_series, endogenous_states_augmented)...)
        if haskey(get_settings(m), :integrated_series)
            integ_series = union(get_setting(m, :integrated_series), integ_series)
        end
        m <= Setting(:integrated_series, integ_series)
    end
    if get_setting(m, :add_cumulative)
        push!(endogenous_states_augmented, setdiff([:y_f_t1, :c_f_t1, :i_f_t1], endogenous_states_augmented)...)
        integ_series = [:cum_z_t, :cum_y_t, :cum_e_gdp_t, :cum_π_t,
                        :cum_y_f_t, :cum_c_t, :cum_c_f_t,
                        :cum_i_t, :cum_i_f_t]
        push!(endogenous_states_augmented, setdiff(integ_series, endogenous_states_augmented)...)
        if haskey(get_settings(m), :integrated_series)
            integ_series = union(get_setting(m, :integrated_series), integ_series)
        end
        m <= Setting(:integrated_series, integ_series)
    end
    if get_setting(m, :add_flexible_price_growth)
        push!(endogenous_states_augmented, setdiff([:y_f_t1, :c_f_t1, :i_f_t1], endogenous_states_augmented)...)
    end
    if haskey(get_settings(m), :add_pgap) ? get_setting(m, :add_pgap) : false
        push!(endogenous_states, setdiff([:pgap_t], endogenous_states)...)
        push!(equilibrium_conditions, setdiff([:eq_pgap], equilibrium_conditions)...)
    end
    if haskey(get_settings(m), :add_ygap) ? get_setting(m, :add_ygap) : false
        push!(endogenous_states, setdiff([:ygap_t], endogenous_states)...)
        push!(equilibrium_conditions, setdiff([:eq_ygap], equilibrium_conditions)...)
    end
    if haskey(get_settings(m), :add_altpolicy_pgap) ? get_setting(m, :add_altpolicy_pgap) : false
        push!(endogenous_states, setdiff([:pgap_t], endogenous_states)...)
        push!(equilibrium_conditions, setdiff([:eq_pgap], equilibrium_conditions)...)
    end
    if haskey(get_settings(m), :add_altpolicy_ygap) ? get_setting(m, :add_altpolicy_ygap) : false
        push!(endogenous_states, setdiff([:ygap_t], endogenous_states)...)
        push!(equilibrium_conditions, setdiff([:eq_ygap], equilibrium_conditions)...)
    end
    if (haskey(get_settings(m), :add_initialize_pgap_ygap_pseudoobs) ? get_setting(m, :add_initialize_pgap_ygap_pseudoobs) : false)
        push!(exogenous_shocks, setdiff([:pgap_sh, :ygap_sh], exogenous_shocks)...)
    end
    if haskey(get_settings(m), :add_rw) ? get_setting(m, :add_rw) : false
        push!(endogenous_states, setdiff([:rw_t, :Rref_t], endogenous_states)...)
        push!(equilibrium_conditions, setdiff([:eq_rw, :eq_Rref], equilibrium_conditions)...)
    end

    if haskey(get_settings(m), :add_iid_cond_obs_gdp_meas_err) ?
        get_setting(m, :add_iid_cond_obs_gdp_meas_err) : false
        push!(endogenous_states_augmented, :e_condgdp_t)
        push!(exogenous_shocks, :condgdp_sh)
    end

    if haskey(get_settings(m), :add_iid_anticipated_obs_gdp_meas_err) ?
        get_setting(m, :add_iid_anticipated_obs_gdp_meas_err) : false
        push!(endogenous_states_augmented, :e_gdpexp_t)
        push!(exogenous_shocks, :gdpexp_sh)
    end

    if haskey(get_settings(m), :add_iid_cond_obs_corepce_meas_err) ?
        get_setting(m, :add_iid_cond_obs_corepce_meas_err) : false
        push!(endogenous_states_augmented, :e_condcorepce_t)
        push!(exogenous_shocks, :condcorepce_sh)
    end

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

function Model1002(subspec::String = "ss10";
                   custom_settings::Dict{Symbol, Setting} = Dict{Symbol, Setting}(),
                   testing = false)

    # Model-specific specifications
    spec               = split(basename(@__FILE__),'.')[1]
    subspec            = subspec
    settings           = Dict{Symbol,Setting}()
    test_settings      = Dict{Symbol,Setting}()
    rng                = MersenneTwister(0)

    # Initialize empty model
    m = Model1002{Float64}(
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
init_parameters!(m::Model1002)
```

Initializes the model's parameters, as well as empty values for the steady-state
parameters (in preparation for `steadystate!(m)` being called to initialize
those).
"""
function init_parameters!(m::Model1002)
    m <= parameter(:α, 0.1596, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Normal(0.30, 0.05), fixed=false,
                   description="α: Capital elasticity in the intermediate goods sector's production function (also known as the capital share).",
                   tex_label="\\alpha")

    m <= parameter(:ζ_p, 0.8940, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.1), fixed=false,
                   description="ζ_p: The Calvo parameter. In every period, intermediate goods producers optimize prices with probability (1-ζ_p). With probability ζ_p, prices are adjusted according to a weighted average of the previous period's inflation (π_t1) and steady-state inflation (π_star).",
                   tex_label="\\zeta_p")

    m <= parameter(:ι_p, 0.1865, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.15), fixed=false,
                   description="ι_p: The weight attributed to last period's inflation in price indexation. (1-ι_p) is the weight attri
buted to steady-state inflation.",
                   tex_label="\\iota_p")

    m <= parameter(:δ, 0.025, fixed=true,
                   description="δ: The capital depreciation rate.",
                   tex_label="\\delta" )

    m <= parameter(:Upsilon, 1.000, (0., 10.), (1e-5, 0.), ModelConstructors.Exponential(), GammaAlt(1., 0.5), fixed=true,
                   description="Υ: The trend evolution of the price of investment goods relative to consumption goods. Set equal to 1.",
                   tex_label="\\Upsilon")

    m <= parameter(:Φ, 1.1066, (1., 10.), (1.00, 10.00), ModelConstructors.Exponential(), Normal(1.25, 0.12), fixed=false,
                   description="Φ: Fixed costs.",
                   tex_label="\\Phi_p")

    m <= parameter(:S′′, 2.7314, (-15., 15.), (-15., 15.), ModelConstructors.Untransformed(), Normal(4., 1.5), fixed=false,
                   description="S'': The second derivative of households' cost of adjusting investment.",
                   tex_label="S''")

    m <= parameter(:h, 0.5347, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.7, 0.1), fixed=false,
                   description="h: Consumption habit persistence.",
                   tex_label="h")

    m <= parameter(:ppsi, 0.6862, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.15), fixed=false,
                   description="ppsi: Utilization costs.",
                   tex_label="\\psi")

    m <= parameter(:ν_l, 2.5975, (1e-5, 10.), (1e-5, 10.), ModelConstructors.Exponential(), Normal(2, 0.75), fixed=false,
                   description="ν_l: The coefficient of relative risk aversion on the labor term of households' utility function.",
                   tex_label="\\nu_l")

    m <= parameter(:ζ_w, 0.9291, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.1), fixed=false,
                   description="ζ_w: (1-ζ_w) is the probability with which households can freely choose wages in each period. With probability ζ_w, wages increase at a geometrically weighted average of the steady state rate of wage increases and last period's productivity times last period's inflation.",
                   tex_label="\\zeta_w")

    m <= parameter(:ι_w, 0.2992, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.15), fixed=false,
                   description="ι_w: The weight attributed to last period's wage in wage indexation. (1-ι_w) is the weight attributed to steady-state wages.",
                   tex_label="\\iota_w")

    m <= parameter(:λ_w, 1.5000, fixed=true,
                   description="λ_w: The wage markup, which affects the elasticity of substitution between differentiated labor services.",
                   tex_label="\\lambda_w")

    m <= parameter(:β, 0.1402, (1e-5, 10.), (1e-5, 10.), ModelConstructors.Exponential(), GammaAlt(0.25, 0.1), fixed=false,
                   scaling = x -> 1/(1 + x/100),
                   description="β: Discount rate.",
                   tex_label="100(\\beta^{-1} - 1)")

    m <= parameter(:ψ1, 1.3679, (1e-5, 10.), (1e-5, 10.00), ModelConstructors.Exponential(), Normal(1.5, 0.25), fixed=false,
                   description="ψ₁: Weight on inflation gap in monetary policy rule.",
                   tex_label="\\psi_1")

    m <= parameter(:ψ2, 0.0388, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.05), fixed=false,
                   description="ψ₂: Weight on output gap in monetary policy rule.",
                   tex_label="\\psi_2")

    m <= parameter(:ψ3, 0.2464, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.05), fixed=false,
                   description="ψ₃: Weight on rate of change of output gap in the monetary policy rule.",
                   tex_label="\\psi_3")

    m <= parameter(:π_star, 0.5000, (1e-5, 10.), (1e-5, 10.), ModelConstructors.Exponential(), GammaAlt(0.75, 0.4), fixed=true,
                   scaling = x -> 1 + x/100,
                   description="π_star: The steady-state rate of inflation.",
                   tex_label="\\pi_*")

    m <= parameter(:σ_c, 0.8719, (1e-5, 10.), (1e-5, 10.), ModelConstructors.Exponential(), Normal(1.5, 0.37), fixed=false,
                   description="σ_c: Coefficient of relative risk aversion.",
                   tex_label="\\sigma_{c}")

    m <= parameter(:ρ, 0.7126, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.75, 0.10), fixed=false,
                   description="ρ: The degree of inertia in the monetary policy rule.",
                   tex_label="\\rho_R")

    m <= parameter(:ϵ_p, 10.000, fixed=true,
                   description="ϵ_p: Curvature parameter in the Kimball aggregator for prices.",
                   tex_label="\\epsilon_{p}")

    m <= parameter(:ϵ_w, 10.000, fixed=true,
                   description="ϵ_w: Curvature parameter in the Kimball aggregator for wages.",
                   tex_label="\\epsilon_{w}")


    # financial frictions parameters
    m <= parameter(:Fω, 0.0300, (1e-5, 0.99999), (1e-5, 0.99), ModelConstructors.SquareRoot(), BetaAlt(0.03, 0.01), fixed=true,
                   scaling = x -> 1 - (1-x)^0.25,
                   description="F(ω): The cumulative distribution function of ω (idiosyncratic iid shock that increases or decreases entrepreneurs' capital).",
                   tex_label="F(\\bar{\\omega})")

    m <= parameter(:spr, 1.7444, (0., 100.), (1e-5, 0.), ModelConstructors.Exponential(), GammaAlt(2., 0.1), fixed=false,
                   scaling = x -> (1 + x/100)^0.25,
                   description="spr_*: Steady-state level of spread.",
                   tex_label="SP_*")

    m <= parameter(:ζ_spb, 0.0559, (1e-5, 0.99999), (1e-5, 0.99), ModelConstructors.SquareRoot(), BetaAlt(0.05, 0.005), fixed=false,
                   description="ζ_spb: The elasticity of the expected exess return on capital (or 'spread') with respect to leverage.",
                   tex_label="\\zeta_{sp,b}")

    m <= parameter(:γ_star, 0.9900, (1e-5, 0.99999), (1e-5, 0.99), ModelConstructors.SquareRoot(), BetaAlt(0.99, 0.002), fixed=true,
                   description="γ_star: Fraction of entrepreneurs who survive and continue operating for another period.",
                   tex_label="\\gamma_*")

    # exogenous processes - level
    m <= parameter(:γ, 0.3673, (-5.0, 5.0), (-5., 5.), ModelConstructors.Untransformed(), Normal(0.4, 0.1), fixed=false,
                   scaling = x -> x/100,
                   description="γ: The log of the steady-state growth rate of technology.",
                   tex_label="100\\gamma")

    m <= parameter(:Lmean, -45.9364, (-1000., 1000.), (-1e3, 1e3), ModelConstructors.Untransformed(), Normal(-45., 5.), fixed=false,
                   description="Lmean: Mean level of hours.",
                   tex_label="\\bar{L}")

    m <= parameter(:g_star, 0.1800, fixed=true,
                   description="g_star: 1 - (c_star + i_star)/y_star.",
                   tex_label="g_*")

    # exogenous processes - autocorrelation
    m <= parameter(:ρ_g, 0.9863, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed=false,
                   description="ρ_g: AR(1) coefficient in the government spending process.",
                   tex_label="\\rho_g")

    m <= parameter(:ρ_b, 0.9410, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed=false,
                   description="ρ_b: AR(1) coefficient in the intertemporal preference shifter process.",
                   tex_label="\\rho_b")

    m <= parameter(:ρ_μ, 0.8735, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed=false,
                   description="ρ_μ: AR(1) coefficient in capital adjustment cost process.",
                   tex_label="\\rho_{\\mu}")

    m <= parameter(:ρ_ztil, 0.9446, (0., 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed=false,
                   description="ρ_ztil: AR(1) coefficient in the technology process.",
                   tex_label="\\rho_{\\tilde{z}}")

    m <= parameter(:ρ_λ_f, 0.8827, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed=false,
                   description="ρ_λ_f: AR(1) coefficient in the price mark-up shock process.",
                   tex_label="\\rho_{\\lambda_f}")

    m <= parameter(:ρ_λ_w, 0.3884, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed=false,
                   description="ρ_λ_w: AR(1) coefficient in the wage mark-up shock process.",
                   tex_label="\\rho_{\\lambda_w}")

    # monetary policy shock - see eqcond
    m <= parameter(:ρ_rm, 0.2135, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed=false,
                   description="ρ_rm: AR(1) coefficient in the monetary policy shock process.",
                   tex_label="\\rho_{r^m}")

    m <= parameter(:ρ_σ_w, 0.9898, (1e-5, 0.99999), (1e-5, 0.99999), ModelConstructors.SquareRoot(), BetaAlt(0.75, 0.15), fixed=false,
                   description="ρ_σ_w: The standard deviation of entrepreneurs' capital productivity follows an exogenous process with mean ρ_σ_w. Innovations to the process are called _spread shocks_.",
                   tex_label="\\rho_{\\sigma_\\omega}")

    m <= parameter(:ρ_μ_e, 0.7500, (1e-5, 0.99999), (1e-5, 0.99), ModelConstructors.SquareRoot(), BetaAlt(0.75, 0.15), fixed=true,
                   description="ρ_μ_e: AR(1) coefficient in the exogenous bankruptcy cost process.",
                   tex_label="\\rho_{\\mu_e}")

    m <= parameter(:ρ_γ, 0.7500, (1e-5, 0.99999), (1e-5, 0.99), ModelConstructors.SquareRoot(), BetaAlt(0.75, 0.15), fixed=true,
                   description="ρ_γ: AR(1) coefficient in the process describing the fraction of entrepreneurs surviving period t.",
                   tex_label="\\rho_{\\gamma}")

    m <= parameter(:ρ_π_star, 0.9900, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed=true,
                   description="ρ_π_star: AR(1) coefficient in the time-varying inflation target process.",
                   tex_label="\\rho_{\\pi_*}")

    m <= parameter(:ρ_lr, 0.6936, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed=false,
                   tex_label="\\rho_{10y}")

    m <= parameter(:ρ_z_p, 0.8910, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed=false,
                   description="ρ_z_p: AR(1) coefficient in the process describing the permanent component of productivity.",
                   tex_label="\\rho_{z^p}")

    m <= parameter(:ρ_tfp, 0.1953, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed=false,
                   tex_label="\\rho_{tfp}")

    m <= parameter(:ρ_gdpdef, 0.5379, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed=false,
                   tex_label="\\rho_{gdpdef}")

    m <= parameter(:ρ_corepce, 0.2320, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed=false,
                   tex_label="\\rho_{pce}")

    m <= parameter(:ρ_gdp, 0., (-0.999, 0.999), (-0.999, 0.999), ModelConstructors.SquareRoot(), Normal(0.0, 0.2), fixed=false,
                   tex_label="\\rho_{gdp}")

    m <= parameter(:ρ_gdi, 0., (-0.999, 0.999), (-0.999, 0.999), ModelConstructors.SquareRoot(), Normal(0.0, 0.2), fixed=false,
                   tex_label="\\rho_{gdi}")

    m <= parameter(:ρ_gdpvar, 0., (-0.999, 0.999), (-0.999, 0.999), ModelConstructors.SquareRoot(), Normal(0.0, 0.4), fixed=false,
                   tex_label="\\varrho_{gdp}")

    m <= parameter(:me_level, 1., (-0.999, 0.999), (-0.999, 0.999), ModelConstructors.Untransformed(), Normal(0.0, 0.4), fixed=true,
                   description="me_level: Indicator of cointegration of GDP and GDI.",
                   tex_label="\\mathcal{C}_{me}")

    # exogenous processes - standard deviation
    m <= parameter(:σ_g, 2.5230, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_g: The standard deviation of the government spending process.",
                   tex_label="\\sigma_{g}")

    m <= parameter(:σ_b, 0.0292, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_b: The standard deviation of the intertemporal preference shifter process.",
                   tex_label="\\sigma_{b}")

    m <= parameter(:σ_μ, 0.4559, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_μ: The standard deviation of the exogenous marginal efficiency of investment shock process.",
                   tex_label="\\sigma_{\\mu}")

    m <= parameter(:σ_ztil, 0.6742, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_ztil: The standard deviation of the process describing the stationary component of productivity.",
                   tex_label="\\sigma_{\\tilde{z}}")

    m <= parameter(:σ_λ_f, 0.1314, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_λ_f: The mean of the process that generates the price elasticity of the composite good. Specifically, the elasticity is (1+λ_{f,t})/(λ_{f_t}).",
                   tex_label="\\sigma_{\\lambda_f}")

    m <= parameter(:σ_λ_w, 0.3864, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{\\lambda_w}")

    m <= parameter(:σ_r_m, 0.2380, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_r_m: The standard deviation of the monetary policy shock.",
                   tex_label="\\sigma_{r^m}")

    m <= parameter(:σ_σ_ω, 0.0428, (1e-7,100.), (1e-5, 0.), ModelConstructors.Exponential(), RootInverseGamma(4, 0.05), fixed=false,
                   description="σ_σ_ω: The standard deviation of entrepreneurs' capital productivity follows an exogenous process with standard deviation σ_σ_ω.",
                   tex_label="\\sigma_{\\sigma_\\omega}")

    m <= parameter(:σ_μ_e, 0.0000, (0., 100.), (1e-5, 0.), ModelConstructors.Exponential(), RootInverseGamma(4, 0.05), fixed=true,
                   description="σ_μ_e: Exogenous bankrupcy costs follow an exogenous process with standard deviation σ_μ_e.",
                   tex_label="\\sigma_{\\mu_e}")

    m <= parameter(:σ_γ, 0.0000, (0., 100.), (1e-5, 0.), ModelConstructors.Exponential(), RootInverseGamma(4, 0.01), fixed=true,
                   description="σ_γ: The fraction of entrepreneurs surviving period t follows an exogenous process with standard deviation σ_γ.",
                   tex_label="\\sigma_{\\gamma}")

    m <= parameter(:σ_π_star, 0.0269, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(6, 0.03), fixed=false,
                   description="σ_π_star: The standard deviation of the inflation target.",
                   tex_label="\\sigma_{\\pi_*}")

    m <= parameter(:σ_lr, 0.1766, (1e-8,10.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.75), fixed=false,
                   tex_label="\\sigma_{10y}")

    m <= parameter(:σ_z_p, 0.1662, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_z_p: The standard deviation of the shock to the permanent component of productivity.",
                   tex_label="\\sigma_{z^p}")

    m <= parameter(:σ_tfp, 0.9391, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{tfp}")

    m <= parameter(:σ_gdpdef, 0.1575, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{gdpdef}")

    m <= parameter(:σ_corepce, 0.0999, (1e-8, 5.),(1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{pce}")

    m <= parameter(:σ_gdp, 0.1, (1e-8, 5.),(1e-8, 5.),ModelConstructors.Exponential(),RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{gdp}")

    m <= parameter(:σ_gdi, 0.1, (1e-8, 5.),(1e-8, 5.),ModelConstructors.Exponential(),RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{gdi}")

    if subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63"]
        m <= parameter(:ρ_ziid, 0., (0., 0.999), (0., 0.999), ModelConstructors.Untransformed(), BetaAlt(0.5, 0.2), fixed=true,
                       description="ρ_ziid: AR(1) coefficient in the iid component of the technology process.",
                       tex_label="\\rho_{z, iid}")
        m <= parameter(:σ_ziid, 0., (0., 1e2), (0., 1e2), ModelConstructors.Exponential(),
                       RootInverseGamma(2. * (5.)^2 ./ 5., sqrt((5.)^2 + .1)), fixed=false,
                       description="σ_ziid: The standard deviation of the process describing the iid component of productivity.",
                       tex_label="\\sigma_{z, iid}")
        m <= parameter(:ρ_biidc, 0., (0., 0.999), (0., 0.999), ModelConstructors.Untransformed(), BetaAlt(0.5, 0.2), fixed=true,
                       description="ρ_biidc: AR(1) coefficient in the iid component of the preference process.",
                       tex_label="\\rho_{b, iid, c}")
        m <= parameter(:σ_biidc, 0., (0., 1e2), (0., 1e2), ModelConstructors.Exponential(),
                       RootInverseGamma(2. * (4.)^2 ./ 4., sqrt((4.)^2  + .1)), fixed=false, # If σ_φ ∼ RootInverseGamma(ν, τ), then σ_φ² ∼ InverseGamma(ν/2, ντ²/2), with mode M given by ν (τ² - M²) = 2 * M²
                       description="σ_biidc: The standard deviation of the process describing the iid component of preferences.",
                       tex_label="\\sigma_{b, iid, c}")
        m <= parameter(:ρ_φ, 0., (0., 0.999), (0., 0.999), ModelConstructors.Untransformed(), BetaAlt(0.5, 0.2), fixed=true,
                       description="ρ_φ: AR(1) coefficient in the labor supply preference process.",
                       tex_label="\\rho_{\\varphi}")
        m <= parameter(:σ_φ, 0., (0., 1e3), (0., 1e3), ModelConstructors.Exponential(),
                       RootInverseGamma(2 * (400.0)^2 ./ (400. * 10), sqrt(1. + (400.0)^2)), fixed=false,
                       description="σ_φ: The standard deviation of the process describing the labor supply preference.",
                       tex_label="\\sigma_{\\varphi}") # If σ_φ ∼ RootInverseGamma(ν, τ), then σ_φ² ∼ InverseGamma(ν/2, ντ²/2)
    end

    if subspec(m) in ["ss62"] # TODO: add ss63
        m <= parameter(:damp_standard_shocks, 1., (0., 1e3), (0., 1e3), fixed=false, # will be fixed later in subspec!(m)
                       description="damp_standard_shocks: Damping factor for standard business cycle shocks during COVID-19.",
                       tex_label="damp standard shocks")
        m <= parameter(:amplify_σ_r_m, 1., (0., 1e3), (0., 1e3), fixed=false, # will be fixed later in subspec!(m)
                       description="amplify_σ_r_m: Amplification factor for monetary policy shock during COVID-19.",
                       tex_label="amplify \\sigma_{r^m}")
        m <= parameter(:amplify_inflation_me, 1., (0., 1e3), (0., 1e3), fixed=false, # will be fixed later in subspec!(m)
                       description="amplify_inflation_me: Amplification factor for inflation measurement error during COVID-19.",
                       tex_label="amplify inflation measurement error")
    end

    if haskey(get_settings(m), :add_initialize_pgap_ygap_pseudoobs) ?
        get_setting(m, :add_initialize_pgap_ygap_pseudoobs) : false
        m <= parameter(:σ_pgap, 0., (0., 1e2), (0., 5.), ModelConstructors.Exponential(),
                       RootInverseGamma(2. * (20.)^2 ./ .1, sqrt((4.)^2  + .1)), fixed=true,
                       tex_label="\\sigma_{pgap}")
        m <= parameter(:σ_ygap, 0., (0., 1e2), (0., 5.), ModelConstructors.Exponential(),
                       RootInverseGamma(2. * (20.)^2 ./ .1, sqrt((4.)^2  + .1)), fixed=true,
                       tex_label="\\sigma_{ygap}")
    end

    if haskey(get_settings(m), :add_iid_cond_obs_gdp_meas_err) ?
        get_setting(m, :add_iid_cond_obs_gdp_meas_err) : false
        m <= parameter(:ρ_condgdp, 0., (-0.999, 0.999), (-0.999, 0.999),
                       ModelConstructors.SquareRoot(), Normal(0.0, 0.2), fixed=true,
                       tex_label="\\rho_{cond gdp}")
        m <= parameter(:σ_condgdp, 0.1, (0., 5.), (0., 5.), ModelConstructors.Exponential(),
                       RootInverseGamma(2, 0.10), fixed=false,
                       tex_label="\\sigma_{cond gdp}")
    end

    if haskey(get_settings(m), :add_iid_anticipated_obs_gdp_meas_err) ?
        get_setting(m, :add_iid_anticipated_obs_gdp_meas_err) : false
        m <= parameter(:ρ_gdpexp, 0., (-0.999, 0.999), (-0.999, 0.999),
                       ModelConstructors.SquareRoot(), Normal(0.0, 0.2), fixed=true,
                       tex_label="\\rho_{gdpexp}")
        m <= parameter(:σ_gdpexp, 0.1, (0., 5.), (0., 5.), ModelConstructors.Exponential(),
                       RootInverseGamma(2, 0.10), fixed=false,
                       tex_label="\\sigma_{gdpexp}")
    end

    if haskey(get_settings(m), :add_iid_cond_obs_corepce_meas_err) ?
        get_setting(m, :add_iid_cond_obs_corepce_meas_err) : false
        m <= parameter(:ρ_condcorepce, 0., (-0.999, 0.999), (-0.999, 0.999),
                       ModelConstructors.SquareRoot(), Normal(0.0, 0.2), fixed=true,
                       tex_label="\\rho_{cond corepce}")
        m <= parameter(:σ_condcorepce, 0.0999, (0., 5.), (0., 5.), ModelConstructors.Exponential(),
                       RootInverseGamma(2, 0.10), fixed=false,
                       tex_label="\\sigma_{cond corepce}")
    end

    for key in get_setting(m, :proportional_antshocks)
        propant_tex_label = DSGE.detexify(key) == key ? string(key) : "\\" * string(DSGE.detexify(key))
        m <= parameter(Symbol(:σ_, key, :_prop), 0., (0., 1e3), (0., 0.), ModelConstructors.Exponential(),
                       RootInverseGamma(2, 1.), fixed=false,
                       description="σ_$(key)_prop: proportional of anticipated shock to contemporaneous shock to $key",
                       tex_label="\\sigma_{$(propant_tex_label)}^{prop}")
    end

    # standard deviations of the anticipated policy shocks
    for i = 1:n_mon_anticipated_shocks_padding(m)
        if i <= n_mon_anticipated_shocks(m)
            m <= parameter(Symbol("σ_r_m$i"), .2, (1e-7, 100.), (1e-5, 0.), ModelConstructors.Exponential(),
                           RootInverseGamma(4, .2), fixed=false,
                           description="σ_r_m$i: Standard deviation of the $i-period-ahead anticipated policy shock.",
                           tex_label=@sprintf("\\sigma_{ant%d}",i))
        else
            m <= parameter(Symbol("σ_r_m$i"), .0, (1e-7, 100.), (1e-5, 0.), ModelConstructors.Exponential(), RootInverseGamma(4, .2), fixed=true,
                           description="σ_r_m$i: Standard deviation of the $i-period-ahead anticipated policy shock.",
                           tex_label=@sprintf("\\sigma_{ant%d}",i))
        end
    end

    for (sh, ant_num) in get_setting(m, :antshocks)
        ant_tex_label = DSGE.detexify(sh) == sh ? string(sh) : "\\" * string(DSGE.detexify(sh))
        for i in 1:ant_num
            m <= parameter(Symbol("σ_$(sh)$i"), 0., (0., 1e3), (0., 1e3), ModelConstructors.Exponential(), RootInverseGamma(4, .2),
                           fixed=false,
                           description="σ_$(sh)$i: Standard deviation of the $i-period-ahead anticipated policy shock.",
                           tex_label="\\sigma_{$(ant_tex_label), ant$(i)}")
        end
    end

    m <= parameter(:η_gz, 0.8400, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.50, 0.20), fixed=false,
                   description="η_gz: Correlate g and z shocks.",
                   tex_label="\\eta_{gz}")

    m <= parameter(:η_λ_f, 0.7892, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.50, 0.20), fixed=false,
                   description="η_λ_f: Moving average component in the price markup shock.",
                   tex_label="\\eta_{\\lambda_f}")

    m <= parameter(:η_λ_w, 0.4226, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.50, 0.20), fixed=false,
                   description="η_λ_w: Moving average component in the wage markup shock.",
                   tex_label="\\eta_{\\lambda_w}")

    m <= parameter(:Iendoα, 0.0000, (0.000, 1.000), (0., 0.), ModelConstructors.Untransformed(), BetaAlt(0.50, 0.20), fixed=true,
                   description="Iendoα: Indicates whether to use the model's endogenous α in the capacity utilization adjustment of total factor productivity.",
                   tex_label="I\\{\\alpha^{model}\\}")

    m <= parameter(:Γ_gdpdef, 1.0354, (-10., 10.), (-10., -10.), ModelConstructors.Untransformed(), Normal(1.00, 2.), fixed=false,
                   tex_label="\\gamma_{gdpdef}")

    m <= parameter(:δ_gdpdef, 0.0181, (-10., 10.), (-10., -10.), ModelConstructors.Untransformed(), Normal(0.00, 2.), fixed=false,
                   tex_label="\\delta_{gdpdef}")

    m <= parameter(:γ_gdi, 1., (-10., 10.), (-10., -10.), ModelConstructors.Untransformed(), Normal(1., 2.), fixed=false,
                   tex_label="\\gamma_{gdi}")

    m <= parameter(:δ_gdi, 0., (-10., 10.), (-10., -10.), ModelConstructors.Untransformed(), Normal(0.00, 2.), fixed=false,
                   tex_label="\\delta_{gdi}")

    # steady states
    m <= SteadyStateParameter(:z_star, NaN, tex_label="\\z_*")
    m <= SteadyStateParameter(:rstar, NaN, tex_label="\\r_*")
    m <= SteadyStateParameter(:Rstarn, NaN, tex_label="\\R_*_n")
    m <= SteadyStateParameter(:r_k_star, NaN, description="Steady-state short-term rate of return on capital.", tex_label="\\r^k_*")
    m <= SteadyStateParameter(:wstar, NaN, tex_label="\\w_*")
    m <= SteadyStateParameter(:Lstar, NaN, tex_label="\\L_*")
    m <= SteadyStateParameter(:kstar, NaN, description="Effective capital that households rent to firms in the steady state.", tex_label="\\k_*")
    m <= SteadyStateParameter(:kbarstar, NaN, description="Total capital owned by households in the steady state.", tex_label="\\bar{k}_*")
    m <= SteadyStateParameter(:istar, NaN, description="Detrended steady-state investment", tex_label="\\i_*")
    m <= SteadyStateParameter(:ystar, NaN, tex_label="\\y_*")
    m <= SteadyStateParameter(:cstar, NaN, tex_label="\\c_*")
    m <= SteadyStateParameter(:wl_c, NaN, tex_label="\\wl_c")
    if subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63"]
        m <= SteadyStateParameter(:φstar, NaN, tex_label="\\varphi_*")
    end
    m <= SteadyStateParameter(:nstar, NaN, tex_label="\\n_*")
    m <= SteadyStateParameter(:vstar, NaN, tex_label="\\v_*")
    m <= SteadyStateParameter(:ζ_spσ_ω, NaN, tex_label="\\zeta_{sp_{\\sigma_\\omega}}")
    m <= SteadyStateParameter(:ζ_spμ_e, NaN, tex_label="\\zeta_{sp_{\\mu_e}}")
    m <= SteadyStateParameter(:ζ_nRk, NaN, tex_label="\\zeta_{n_R_k}")
    m <= SteadyStateParameter(:ζ_nR, NaN, tex_label="\\zeta_{n_R}")
    m <= SteadyStateParameter(:ζ_nqk, NaN, tex_label="\\zeta_{n_q_k}")
    m <= SteadyStateParameter(:ζ_nn, NaN, tex_label="\\zeta_{nn}")
    m <= SteadyStateParameter(:ζ_nμ_e, NaN, tex_label="\\zeta_{n_{\\mu_e}}")
    m <= SteadyStateParameter(:ζ_nσ_ω, NaN, tex_label="\\zeta_{n_{\\sigma_\\omega}}")
end

"""
```
steadystate!(m::Model1002)
```

Calculates the model's steady-state values. `steadystate!(m)` must be called whenever
the parameters of `m` are updated.
"""
function steadystate!(m::Model1002)
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

    if subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63"]
        m[:φstar] = 0. # log(1) in steady state
    end

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
    ωbarstar    = ω_fn(zω_star, σ_ω_star)

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
    if subspec(m) in ["ss2", "ss8"]
        wekstar = (1-m[:γ_star]/m[:β])*nkstar - m[:γ_star]/m[:β]*(m[:spr]*(1-μ_estar*Gstar) - 1)
    else
        if subspec(m) == "ss9"
            betabar_inverse = exp( (σ_ω_star -1) * m[:z_star]) / m[:β]
        else
            betabar_inverse = exp( (m[:σ_c] -1) * m[:z_star]) / m[:β]
        end
        wekstar = (1-(m[:γ_star]*betabar_inverse))*nkstar - m[:γ_star]*betabar_inverse*(m[:spr]*(1-μ_estar*Gstar) - 1)
    end
    vkstar      = (nkstar-wekstar)/m[:γ_star]

    # evaluate nstar and vstar
    if subspec(m) in ["ss2", "ss8"]
        m[:nstar] = nkstar*m[:kstar]
        m[:vstar] = vkstar*m[:kstar]
    else
        m[:nstar] = nkstar*m[:kbarstar]
        m[:vstar] = vkstar*m[:kbarstar]
    end

    # a couple of combinations
    ΓμG         = Γstar - μ_estar*Gstar
    ΓμGprime    = dΓdω_star - μ_estar*dGdω_star

    # elasticities wrt ωbar
    ζ_bw        = ζ_bω_fn(zω_star, σ_ω_star, m[:spr])
    ζ_zw        = ζ_zω_fn(zω_star, σ_ω_star, m[:spr])
    ζ_bw_zw     = ζ_bw/ζ_zw

    # elasticities wrt σ_ω
    ζ_bσ_ω      = σ_ω_star * (((1 - μ_estar*dGdσstar/dΓdσstar) /
        (1 - μ_estar*dGdω_star/dΓdω_star) - 1)*dΓdσstar*m[:spr] + μ_estar*nkstar*
            (dGdω_star*d2Γdωdσstar - dΓdω_star*d2Gdωdσstar)/ΓμGprime^2) /
                ((1 - Γstar)*m[:spr] + dΓdω_star/ΓμGprime*(1-nkstar))
    ζ_zσ_ω      = σ_ω_star * (dΓdσstar - μ_estar*dGdσstar) / ΓμG
    m[:ζ_spσ_ω] = (ζ_bw_zw*ζ_zσ_ω - ζ_bσ_ω) / (1-ζ_bw_zw)

    # elasticities wrt μ_e
    if subspec(m) in ["ss2", "ss8"]
        ζ_bμ_e  = μ_estar * (nkstar*dΓdω_star*dGdω_star/ΓμGprime+dΓdω_star*Gstar*m[:spr]) /
            ((1-Γstar)*ΓμGprime*m[:spr] + dΓdω_star*(1-nkstar))
    else
        ζ_bμ_e  = -μ_estar * (nkstar*dΓdω_star*dGdω_star/ΓμGprime+dΓdω_star*Gstar*m[:spr]) /
            ((1-Γstar)*ΓμGprime*m[:spr] + dΓdω_star*(1-nkstar))
    end
    ζ_zμ_e      = -μ_estar*Gstar/ΓμG
    m[:ζ_spμ_e] = (ζ_bw_zw*ζ_zμ_e - ζ_bμ_e) / (1-ζ_bw_zw)

    # some ratios/elasticities
    Rkstar      = m[:spr]*m[:π_star]*m[:rstar] # (r_k_star+1-δ)/Upsilon*π_star
    ζ_gw        = dGdω_star/Gstar*ωbarstar
    ζ_Gσ_ω      = dGdσstar/Gstar*σ_ω_star

    # elasticities for the net worth evolution
    m[:ζ_nRk]     = m[:γ_star]*Rkstar/m[:π_star]/exp(m[:z_star])*(1+Rhostar)*(1 - μ_estar*Gstar*(1 - ζ_gw/ζ_zw))
    if subspec(m) in ["ss2", "ss8"]
        m[:ζ_nR]  = m[:γ_star]/m[:β]*(1+Rhostar)*(1 - nkstar + μ_estar*Gstar*m[:spr]*ζ_gw/ζ_zw)
        m[:ζ_nqk] = m[:γ_star]*Rkstar/m[:π_star]/exp(m[:z_star])*(1+Rhostar)*(1 - μ_estar*Gstar*(1+ζ_gw/ζ_zw/Rhostar)) - m[:γ_star]/m[:β]*(1+Rhostar)
        m[:ζ_nn]  = m[:γ_star]/m[:β] + m[:γ_star]*Rkstar/m[:π_star]/exp(m[:z_star])*(1+Rhostar)*μ_estar*Gstar*ζ_gw/ζ_zw/Rhostar
    else
        m[:ζ_nR]  = m[:γ_star]*betabar_inverse*(1+Rhostar)*(1 - nkstar + μ_estar*Gstar*m[:spr]*ζ_gw/ζ_zw)
        m[:ζ_nqk] = m[:γ_star]*Rkstar/m[:π_star]/exp(m[:z_star])*(1+Rhostar)*(1 - μ_estar*Gstar*(1+ζ_gw/ζ_zw/Rhostar)) - m[:γ_star]*betabar_inverse*(1+Rhostar)
        m[:ζ_nn]  = m[:γ_star]*betabar_inverse + m[:γ_star]*Rkstar/m[:π_star]/exp(m[:z_star])*(1+Rhostar)*μ_estar*Gstar*ζ_gw/ζ_zw/Rhostar
    end

    m[:ζ_nμ_e]    = m[:γ_star]*Rkstar/m[:π_star]/exp(m[:z_star])*(1+Rhostar)*μ_estar*Gstar*(1 - ζ_gw*ζ_zμ_e/ζ_zw)
    m[:ζ_nσ_ω]    = m[:γ_star]*Rkstar/m[:π_star]/exp(m[:z_star])*(1+Rhostar)*μ_estar*Gstar*(ζ_Gσ_ω-ζ_gw/ζ_zw*ζ_zσ_ω)

    return m
end

function model_settings!(m::Model1002)

    default_settings!(m)

    # Anticipated shocks
    m <= Setting(:n_mon_anticipated_shocks, 6,
                 "Number of anticipated policy shocks")
    m <= Setting(:n_mon_anticipated_shocks_padding, 20,
                 "Padding for anticipated policy shocks")

    m <= Setting(:n_z_anticipated_shocks, 0,
                 "Number of anticipated policy shocks")
    m <= Setting(:n_z_anticipated_shocks_padding, 0,
                 "Padding for anticipated policy shocks")

    m <= Setting(:ant_eq_mapping,  Dict{Symbol, Symbol}(:z => :ztil),
                 "Dictionary mapping name of anticipated shock to name of state variable for exogenous process affected by the shock")
    m <= Setting(:ant_eq_E_mapping, Dict{Symbol, Symbol}(),
                 "Dictionary mapping name of anticipated shock to name of state variable for" *
                 " one-period ahead expectation of exogenous process affected by the shock") # This mapping is intended only for anticipated shocks that are not scaled by a constant, i.e. the coefficient in the Γ0 eqcond matrix is -1.0
    m <= Setting(:proportional_antshocks, Vector{Symbol}(undef, 0),
                 "Vector of anticipated shocks that are proportional to the contemporaneous shock")
    m <= Setting(:contemporaneous_and_proportional_antshocks, Vector{Symbol}(undef, 0),
                 "Vector of anticipated shocks which are both contemporaneous and shocks proportional to contemporaneous shocks")
    m <= Setting(:antshocks, Dict{Symbol, Int}(),
                 "Dictionary mapping name of anticipated shock to the number of periods of anticipation.")

    # Data
    m <= Setting(:data_id, 3, "Dataset identifier")
    if get_setting(m, :cond_id) in collect(1:5)
        m <= Setting(:cond_full_names, [:obs_gdp, :obs_corepce, :obs_spread, :obs_nominalrate, :obs_longrate],
                     "Observables used in conditional forecasts")
    elseif get_setting(m, :cond_id) == 6
        m <= Setting(:cond_full_names, [:obs_gdp, :obs_corepce, :obs_spread, :obs_nominalrate, :obs_longrate,
                                        :obs_gdpdeflator], "Observables used in conditional forecasts")
    end
    m <= Setting(:cond_semi_names, [:obs_spread, :obs_nominalrate, :obs_longrate],
                 "Observables used in semiconditional forecasts")
    if subspec(m) in ["ss16", "ss17"]
        m <= Setting(:laborshare_base_period, DSGE.quartertodate("1964-Q1"), "Base year for labor share series to provide an initial condition")
    end

    # Forecast
    m <= Setting(:use_population_forecast, true,
                 "Whether to use population forecasts as data")
    m <= Setting(:shockdec_startdate, Nullable(quartertodate("2007-Q1")),
                 "Date of start of shock decomposition output period. If null, then shockdec starts at date_mainsample_start")

    # COVID-19 settings
    if subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63"]
        m <= Setting(:antshocks, Dict{Symbol, Int}(:biidc => 1, :φ => 1, :ziid => 1))
        m <= Setting(:ant_eq_mapping, Dict{Symbol, Symbol}(:biidc => :biidc, :φ => :φ, :ziid => :ziid))
        m <= Setting(:ant_eq_E_mapping, Dict{Symbol, Symbol}(:φ => :Eφ))
        m <= Setting(:proportional_antshocks, Symbol[:biidc, :φ, :ziid])
    end
    if subspec(m) in ["ss62", "ss63"]
        m <= Setting(:add_pseudo_gdp, true),
        m <= Setting(:add_pseudo_corepce, true)
        m <= Setting(:add_anticipated_obs_gdp, true)
        m <= Setting(:n_anticipated_obs_gdp, 1)
        m <= Setting(:filename_anticipated_obs_gdp, "MEDIANANTGDP")
        m <= Setting(:contemporaneous_and_proportional_antshocks, Symbol[:biidc])
        m <= Setting(:meas_err_anticipated_obs_gdp, 1.)
        m <= Setting(:add_iid_anticipated_obs_gdp_meas_err, false)
        m <= Setting(:add_initialize_pgap_ygap_pseudoobs, true)
    end

    # Add AIT for 2020-Q3 on
    m <= Setting(:flexible_ait_policy_change_date, Date(2020, 9, 30))
    if subspec(m) in ["ss30", "ss59", "ss60", "ss61"]
        m <= Setting(:flexible_ait_policy_change, false,
                     "Indicator for whether 2020-Q3 switch in monetary policy rule to AIT is on")
        m <= Setting(:add_pgap, true)
        m <= Setting(:add_ygap, true)
    elseif subspec(m) in ["ss62", "ss63"]
        m <= Setting(:flexible_ait_policy_change, false,
                     "Indicator for whether 2020-Q3 switch in monetary policy rule to AIT is on")
        m <= Setting(:add_pgap, true)
        m <= Setting(:add_ygap, true)
        m <= Setting(:skip_altpolicy_state_init, true)
    else
        m <= Setting(:flexible_ait_policy_change, false,
                     "Indicator for whether 2020-Q3 switch in monetary policy rule to AIT is on")
        m <= Setting(:add_pgap, false)
        m <= Setting(:add_ygap, false)
    end

    # ZLB via zero-rate
    m <= Setting(:zlb_rule_value, 0.)

    # Additional pseudo-observables and integrated series
    m <= Setting(:add_laborshare_measurement, false)
    m <= Setting(:add_laborproductivity_measurement, false)
    m <= Setting(:add_nominalgdp_level, false)
    m <= Setting(:add_nominalgdp_growth, false)
    m <= Setting(:add_cumulative, false)
    m <= Setting(:add_flexible_price_growth, false)

    nothing
end

"""
```
parameter_groupings(m::Model1002)
```

Returns an `OrderedDict{String, Vector{Parameter}}` mapping descriptions of
parameter groupings (e.g. \"Policy Parameters\") to vectors of
`Parameter`s. This dictionary is passed in as a keyword argument to
`prior_table`.
"""
function parameter_groupings(m::Model1002)
    policy     = [[:ψ1, :ψ2, :ψ3, :ρ, :ρ_rm, :σ_r_m];
                  [Symbol("σ_r_m$i") for i = 1:n_mon_anticipated_shocks(m)]]
    sticky     = [:ζ_p, :ι_p, :ϵ_p, :ζ_w, :ι_w, :ϵ_w]
    other_endo = [:γ, :α, :β, :σ_c, :h, :ν_l, :δ, :Φ, :S′′, :ppsi, :π_star,
                  :Γ_gdpdef, :δ_gdpdef, :Lmean, :λ_w, :g_star]
    financial  = [:Fω, :spr, :ζ_spb, :γ_star]
    processes  = [:ρ_g, :ρ_b, :ρ_μ, :ρ_ztil, :ρ_σ_w, :ρ_π_star, :ρ_z_p, :ρ_λ_f, :ρ_λ_w, :η_λ_f, :η_λ_w,
                  :σ_g, :σ_b, :σ_μ, :σ_ztil, :σ_σ_ω, :σ_π_star, :σ_z_p, :σ_λ_f, :σ_λ_w, :η_gz]
    error      = [:me_level, :ρ_gdp, :ρ_gdi, :ρ_lr, :ρ_tfp, :ρ_gdpdef, :ρ_corepce,
                  :ρ_gdpvar, :σ_gdp, :σ_gdi, :σ_lr, :σ_tfp, :σ_gdpdef, :σ_corepce]

    all_keys     = Vector[policy, sticky, other_endo, financial, processes, error]
    descriptions = ["Policy Parameters", "Nominal Rigidities Parameters",
                    "Other Endogenous Propagation and Steady State Parameters",
                    "Financial Frictions Parameters", "Exogenous Process Parameters",
                    "Measurement Error Parameters"]

    if subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63"]
        covid = [:σ_ziid, :σ_biidc, :σ_φ]
        for (sh, ant_num) in get_setting(m, :antshocks)
            for i in 1:ant_num
                push!(covid, Symbol("σ_$(sh)$i"))
            end
        end
        for key in get_setting(m, :proportional_antshocks)
            push!(covid, Symbol(:σ_, key, :_prop))
        end
        push!(all_keys, covid)
        push!(descriptions, "COVID-19 Parameters")
    end

    all_params = map(keys -> [m[θ]::Parameter for θ in keys], all_keys)
    groupings  = OrderedDict{String, Vector{Parameter}}(zip(descriptions, all_params))

    # Ensure no parameters missing
    incl_params    = vcat(collect(values(groupings))...)
    excl_params_sym = vcat([:Upsilon, :ρ_μ_e, :ρ_γ, :σ_μ_e, :σ_γ, :Iendoα, :γ_gdi, :δ_gdi],
                           [Symbol("σ_r_m$i") for i=n_mon_anticipated_shocks(m)+1:n_mon_anticipated_shocks_padding(m)])
    if subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63"]
        push!(excl_params_sym, :ρ_ziid, :ρ_biidc, :ρ_φ)
    end
    if haskey(get_settings(m), :add_initialize_pgap_ygap_pseudoobs) ?
        get_setting(m, :add_initialize_pgap_ygap_pseudoobs) : false
        push!(excl_params_sym, :σ_pgap, :σ_ygap)
    end
    excl_params = [m[θ] for θ in excl_params_sym]

    @assert isempty(setdiff(m.parameters, vcat(incl_params, excl_params)))

    return groupings
end

"""
```
shock_groupings(m::Model1002)
```

Returns a `Vector{ShockGroup}`, which must be passed in to
`plot_shock_decomposition`. See `?ShockGroup` for details.
"""
function shock_groupings(m::Model1002)
    if subspec(m) in ["ss59", "ss60", "ss61", "ss62", "ss63"]
        gov = ShockGroup("g", [:g_sh], RGB(0.70, 0.13, 0.13)) # firebrick
        bet = ShockGroup("b", [:b_sh], RGB(0.3, 0.3, 1.0))
        fin = ShockGroup("FF", [:γ_sh, :μ_e_sh, :σ_ω_sh], RGB(0.29, 0.0, 0.51)) # indigo
        # tfp = ShockGroup("tfp", [:ztil_sh], RGB(1.0, 0.55, 0.0)) # darkorange
        tfp = ShockGroup("tfp", [:ztil_sh, :zp_sh], RGB(1.0, 0.55, 0.0)) # darkorange
        pmu = ShockGroup("mkp", [:λ_f_sh, :λ_w_sh], RGB(0.60, 0.80, 0.20)) # yellowgreen
        wmu = ShockGroup("w-mkp", [:λ_w_sh], RGB(0.0, 0.5, 0.5)) # teal
        phi = ShockGroup("phi", [:φ_sh], RGB(0.5, 0.5, 0.))
        pol = ShockGroup("pol", vcat([:rm_sh], [Symbol("rm_shl$i") for i = 1:n_mon_anticipated_shocks(m)]),
                         RGB(1.0, 0.84, 0.0)) # gold
        pis = ShockGroup("pi-LR", [:π_star_sh], RGB(1.0, 0.75, 0.793)) # pink
        mei = ShockGroup("mu", [:μ_sh], :cyan)
        if haskey(get_settings(m), :add_iid_cond_obs_gdp_meas_err) ?
            get_setting(m, :add_iid_cond_obs_gdp_meas_err) : false
            mea = ShockGroup("me", [:lr_sh, :tfp_sh, :gdpdef_sh, :corepce_sh, :gdp_sh, :gdi_sh, :condgdp_sh, :condcorepce_sh],
                             RGB(0.0, 0.8, 0.0))
        elseif haskey(get_settings(m), :add_iid_anticipated_obs_gdp_meas_err) ?
            get_setting(m, :add_iid_anticipated_obs_gdp_meas_err) : false
            mea = ShockGroup("me", [:lr_sh, :tfp_sh, :gdpdef_sh, :corepce_sh, :gdp_sh, :gdi_sh, :gdpexp_sh],
                             RGB(0.0, 0.8, 0.0))
        else
            mea = ShockGroup("me", [:lr_sh, :tfp_sh, :gdpdef_sh, :corepce_sh, :gdp_sh, :gdi_sh], RGB(0.0, 0.8, 0.0))
        end
        zpe = ShockGroup("zp", [:zp_sh], RGB(0.0, 0.3, 0.0))
        det = ShockGroup("dt", [:dettrend], :gray40)
        oth = ShockGroup("other", [:dettrend, :g_sh, :π_star_sh, :μ_sh], :gray40) # :dettrend

        # COVID-19 Shocks
        betcovid = ShockGroup("biidc", haskey(m.exogenous_shocks, :biidc_shl1) ? [:biidc_sh, :biidc_shl1] : [:biidc_sh],
                              RGB(0.70, 0.13, 0.13))
        zcovid   = ShockGroup("ziid", [:ziid_sh], RGB(0., 0.5, 0.5))
        φcovid   = ShockGroup("phi", [:φ_sh], RGB(0.5, 0.5, 0.))
        ocovid   = ShockGroup("Other COVID", [:ziid_sh, :φ_sh], RGB(0., 0.5, 0.5))

        if haskey(m.exogenous_shocks, :pgap_sh)
            push!(ocovid.shocks, :pgap_sh)
        end
        if haskey(m.exogenous_shocks, :ygap_sh)
            push!(ocovid.shocks, :ygap_sh)
        end

        # Time-varying CCC, e.g. from temporary ZLB
        if haskey(get_settings(m), :gensys2) ? get_setting(m, :gensys2) : false
            st = ShockGroup("States Trend", [:StatesTrend], :darkgreen) # :dettrend
            return [betcovid, ocovid, bet, fin, tfp, pmu, pol, mea, oth, st]
        else
            return [betcovid, ocovid, bet, fin, tfp, pmu, pol, mea, oth]
        end
        # return [gov, bet, fin, tfp, pmu, wmu, pol, pis, mei, mea, zpe, det]
        # return [betcovid, zcovid, φcovid, bet, fin, tfp, pmu, pol, mea, oth]
    elseif subspec(m) != "ss12"
        gov = ShockGroup("g", [:g_sh], RGB(0.70, 0.13, 0.13)) # firebrick
        bet = ShockGroup("b", [:b_sh], RGB(0.3, 0.3, 1.0))
        fin = ShockGroup("FF", [:γ_sh, :μ_e_sh, :σ_ω_sh], RGB(0.29, 0.0, 0.51)) # indigo
        tfp = ShockGroup("z", [:ztil_sh], RGB(1.0, 0.55, 0.0)) # darkorange
        pmu = ShockGroup("p-mkp", [:λ_f_sh], RGB(0.60, 0.80, 0.20)) # yellowgreen
        wmu = ShockGroup("w-mkp", [:λ_w_sh], RGB(0.0, 0.5, 0.5)) # teal
        pol = ShockGroup("pol", vcat([:rm_sh], [Symbol("rm_shl$i") for i = 1:n_mon_anticipated_shocks(m)]),
                         RGB(1.0, 0.84, 0.0)) # gold
        pis = ShockGroup("pi-LR", [:π_star_sh], RGB(1.0, 0.75, 0.793)) # pink
        mei = ShockGroup("mu", [:μ_sh], :cyan)
        mea = ShockGroup("me", [:lr_sh, :tfp_sh, :gdpdef_sh, :corepce_sh, :gdp_sh, :gdi_sh], RGB(0.0, 0.8, 0.0))
        zpe = ShockGroup("zp", [:zp_sh], RGB(0.0, 0.3, 0.0))
        det = ShockGroup("dt", [:dettrend], :gray40)

        return [gov, bet, fin, tfp, pmu, wmu, pol, pis, mei, mea, zpe, det]
    else
        # financial, productivity, other, measurement errors
        fin = ShockGroup("FF", [:γ_sh, :μ_e_sh, :σ_ω_sh], RGB(0.29, 0.0, 0.51)) # indigo
        prod = ShockGroup("Prod", [:ztil_sh, :zp_sh], RGB(1.0, 0.55, 0.0)) # darkorange
        mea = ShockGroup("Measurement", [:lr_sh, :tfp_sh, :gdpdef_sh, :corepce_sh,
                                         :gdp_sh, :gdi_sh], RGB(0.0, 0.8, 0.0))
        other = ShockGroup("Other", [:g_sh, :b_sh, :λ_f_sh, :λ_w_sh,
                                     vcat([:rm_sh],
                                     [Symbol("rm_shl$i") for i = 1:n_mon_anticipated_shocks(m)])...,
                                     :π_star_sh, :μ_sh,
                                     :dettrend], RGB(0.70, 0.13, 0.13)) # firebrick
        return [fin, prod, mea, other]
    end
end
