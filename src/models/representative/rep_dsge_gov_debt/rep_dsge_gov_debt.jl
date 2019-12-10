# Lightweight Representative Agent Version of HetDSGEGovDebt
mutable struct RepDSGEGovDebt{T} <: AbstractRepModel{T}
    parameters::ParameterVector{T}                # vector of all time-invariant model parameters
    steady_state::ParameterVector{T}              # model steady-state values

    grids::OrderedDict{Symbol,Union{Grid, Array}} # Human-readable names for all the model
    keys::OrderedDict{Symbol,Int}                 # parameters and steady-states

    state_variables::Vector{Symbol}               # Vector of symbols of the state variables
    jump_variables::Vector{Symbol}                # Vector of symbols of the jump variables

    # Vector of unnormalized ranges of indices
    endogenous_states::OrderedDict{Symbol,Int}
    exogenous_shocks::OrderedDict{Symbol,Int}
    expected_shocks::OrderedDict{Symbol,Int}
    equilibrium_conditions::OrderedDict{Symbol,Int}
    endogenous_states_augmented::OrderedDict{Symbol, Int}
    observables::OrderedDict{Symbol,Int}

    spec::String                                  # Model specification number (eg "m990")
    subspec::String                               # Model subspecification (eg "ss0")
    settings::Dict{Symbol,Setting}                # Settings/flags for computation
    testing::Bool                                 # Whether we are in testing mode or not
    observable_mappings::OrderedDict{Symbol, Observable}
end

description(m::RepDSGEGovDebt) = "RepDSGEGovDebt, $(m.subspec)"

function RepDSGEGovDebt(het::HetDSGEGovDebt;
                        custom_settings::Dict{Symbol, Setting} = Dict{Symbol, Setting}())

    @assert !isnan(het[:βstar].value) "Must have already called steadystate! on `het` before instantiating a RepDSGEGovDebt model"

    # Model-specific specifications
    spec               = "rep_dsge"
    subspec            = deepcopy(het.subspec)
    settings           = deepcopy(het.settings)

    # initialize empty model
    m = RepDSGEGovDebt{Float64}(
            # model parameters and steady state values
            Vector{AbstractParameter{Float64}}(),
            Vector{Float64}(),

            # grids and keys
            OrderedDict{Symbol,Union{Grid, Array}}(), OrderedDict{Symbol,Int}(),

            # state_inds, jump_inds
            Vector{Symbol}(), Vector{Symbol}(),

            # model indices
            OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(),
            OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(),
            OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(),

            spec,
            subspec,
            settings,
            deepcopy(het.testing),
            deepcopy(het.observable_mappings))

    # Set settings
    model_settings!(m, het)
    for custom_setting in values(custom_settings)
        m <= custom_setting
    end

    # Initialize model indices
    init_model_indices!(m)

    # Initialize parameters
    init_parameters!(m, het)

    # Initialize grids
    init_grids!(m, het)

    # initialize subspec from het
    init_subspec!(m, het)

    # Solve for the steady state
    steadystate!(m, het, same_beta = get_setting(m, :use_same_beta_as_hank))

    return m
end

function model_settings!(m::RepDSGEGovDebt, het::HetDSGEGovDebt)
    for (key, setting) in deepcopy(het.settings)
        m <= setting
    end

    # Total grid x*s
    # Not sure what to do with the below line. Leave commented for now
    # m <= Setting(:n, (get_setting(m, :nx1) +get_setting(m, :nx2)),
                 # "Total grid size, multiplying across grid dimensions.")
    # The rest of these lines were taken directly from the rank_irfs.jl file
    m <= Setting(:n_jumps, 15, "num scalar jumps")
    m <= Setting(:n_backward_looking_states, 14, "num scalar states")
    m <= Setting(:n_model_states, get_setting(m, :n_backward_looking_states) +
                 get_setting(m, :n_jumps),
                 "Number of 'states' in the state space model. Because backward and forward
                 looking variables need to be explicitly tracked for the Klein solution
                 method, we have n_states and n_jumps")
    m <= Setting(:n_model_states_augmented, get_setting(m, :n_model_states) +
                 length([:c_t1]),
                 "The number of 'states' in the augmented state space model.")

    m <= Setting(:use_same_beta_as_hank, false)
end

function init_model_indices!(m::RepDSGEGovDebt)

    # Endogenous states
    states = collect([:k′_t, :R′_t1, :i′_t1, :y′_t1, :w′_t1, :I′_t1, :bg′_t,
                      :b′_t, :g′_t, :z′_t, :μ′_t, :λ_w′_t, :λ_f′_t, :rm′_t])

    jumps = collect([:R′_t, :i′_t, :t′_t, :w′_t, :L′_t, :π′_t, :π_w′_t, :margutil′_t, :y′_t, :I′_t,
                     :mc′_t, :Q′_t, :capreturn′_t, :l′_t, :tg′_t])

    endogenous_states = collect(vcat(states, jumps))

    # Exogenous shocks
    exogenous_shocks = collect([:b_sh,:g_sh,:z_sh,:μ_sh,:λ_w_sh,:λ_f_sh,:rm_sh])

    # Equilibrium conditions
    equilibrium_conditions = collect([:eq_euler, :eq_market_clearing, :eq_lambda, :eq_transfers,
                                      :eq_investment, :eq_tobin_q, :eq_capital_accumulation,
                                      :eq_wage_phillips, :eq_price_phillips, :eq_marginal_cost,
                                      :eq_gdp, :eq_optimal_kl, :eq_taylor, :eq_fisher,
                                      :eq_nominal_wage_inflation, :eq_fiscal_rule,
                                      :eq_g_budget_constraint, :LR, :LI, :LY,
                                      :LW, :LX, :eq_b, :eq_g, :eq_z, :eq_μ, :eq_λ_w,
                                      :eq_λ_f, :eq_rm]) #, :eq_consumption])

    # Additional states added after solving model
    # Lagged states and observables measurement error
    endogenous_states_augmented = [:c_t1]

    observables = keys(m.observable_mappings)
    endo   = m.endogenous_states
    eqcond = equilibrium_conditions

    for (i,k) in enumerate(endogenous_states);m.endogenous_states[k] = i end
    for (i,k) in enumerate(endogenous_states_augmented); m.endogenous_states_augmented[k] = i+length(endogenous_states) end
    for (i,k) in enumerate(exogenous_shocks); m.exogenous_shocks[k]  = i end
    for (i,k) in enumerate(observables);      m.observables[k]       = i end
    for (i,k) in enumerate(equilibrium_conditions);      m.equilibrium_conditions[k]      = i end
end

function init_parameters!(m::RepDSGEGovDebt, het::HetDSGEGovDebt)
    for param in deepcopy(het.parameters)
        m <= param
    end
end

function init_grids!(m::RepDSGEGovDebt, het::HetDSGEGovDebt)
    m.grids = deepcopy(het.grids)
end

function steadystate!(m::RepDSGEGovDebt)
    ns = get_setting(m, :ns)

    #m <= SteadyStateParameter(:βstar, het[:βstar].value, description = "Steady-state discount factor", tex_label = "\\beta_*")
    #m[:r] = exp(m[:γ])/m[:βstar]
    m <= SteadyStateParameter(:Rkstar, m[:r] + m[:δ],
                              description = "Rental rate on capital", tex_label = "Rk_*")
    m <= SteadyStateParameter(:ωstar,
                              (m[:α]^(m[:α]/(1-m[:α])))*(1-m[:α])*m[:Rkstar]^(-m[:α]/(1-m[:α])),
                              description = "Real wage", tex_label = "\\omega_*")
    m <= SteadyStateParameter(:klstar, (m[:α]/(1-m[:α]))*(m[:ωstar]/m[:Rkstar])*exp(m[:γ]),
                              description = "Capital/Labor ratio", tex_label = "kl_*")
    m <= SteadyStateParameter(:kstar, m[:klstar]*m[:H], description = "Capital",
                              tex_label = "k_*")
    m <= SteadyStateParameter(:ystar, (exp(-m[:α]*m[:γ])*(m[:kstar]^m[:α])*m[:H]^(1-m[:α])),
                              description = "GDP", tex_label = "y_*")
    m <= SteadyStateParameter(:xstar, (1-(1-m[:δ])*exp(-m[:γ]))*m[:kstar],
                              description = "Investment", tex_label = "x_*")
    m <= SteadyStateParameter(:cstar, m[:ystar]/m[:g] - m[:xstar],
                              description = "Steady-state consumption", tex_label = "c_*")
end

function steadystate!(m::RepDSGEGovDebt, het::HetDSGEGovDebt;
                      same_beta::Bool = false)
    ns = get_setting(m, :ns)

    if same_beta
        m <= SteadyStateParameter(:βstar, het[:βstar].value, description = "Steady-state discount factor", tex_label = "\\beta_*")
        m[:r] = exp(m[:γ])/m[:βstar]
        m <= SteadyStateParameter(:Rkstar, m[:r] + m[:δ],
                                  description = "Rental rate on capital", tex_label = "Rk_*")
        m <= SteadyStateParameter(:ωstar,
                                  (m[:α]^(m[:α]/(1-m[:α])))*(1-m[:α])*m[:Rkstar]^(-m[:α]/(1-m[:α])),
                                  description = "Real wage", tex_label = "\\omega_*")
        m <= SteadyStateParameter(:klstar, (m[:α]/(1-m[:α]))*(m[:ωstar]/m[:Rkstar])*exp(m[:γ]),
                                  description = "Capital/Labor ratio", tex_label = "kl_*")
        m <= SteadyStateParameter(:kstar, m[:klstar]*m[:H], description = "Capital",
                                  tex_label = "k_*")
        m <= SteadyStateParameter(:ystar, (exp(-m[:α]*m[:γ])*(m[:kstar]^m[:α])*m[:H]^(1-m[:α])),
                                  description = "GDP", tex_label = "y_*")
        m <= SteadyStateParameter(:xstar, (1-(1-m[:δ])*exp(-m[:γ]))*m[:kstar],
                                  description = "Investment", tex_label = "x_*")
        m <= SteadyStateParameter(:cstar, m[:ystar]/m[:g] - m[:xstar],
                                  description = "Steady-state consumption", tex_label = "c_*")
    else
        m <= SteadyStateParameter(:βstar, exp(m[:γ])/(1 + m[:r]), description = "Steady-state discount factor",
                                  tex_label = "\\beta_*")
        m <= SteadyStateParameter(:Rkstar, m[:r] + m[:δ],
                                  description = "Rental rate on capital", tex_label = "Rk_*")
        m <= SteadyStateParameter(:ωstar,
                                  (m[:α]^(m[:α]/(1-m[:α])))*(1-m[:α])*m[:Rkstar]^(-m[:α]/(1-m[:α])),
                                  description = "Real wage", tex_label = "\\omega_*")
        m <= SteadyStateParameter(:klstar, (m[:α]/(1-m[:α]))*(m[:ωstar]/m[:Rkstar])*exp(m[:γ]),
                                  description = "Capital/Labor ratio", tex_label = "kl_*")
        m <= SteadyStateParameter(:kstar, m[:klstar]*m[:H], description = "Capital",
                                  tex_label = "k_*")
        m <= SteadyStateParameter(:ystar, (exp(-m[:α]*m[:γ])*(m[:kstar]^m[:α])*m[:H]^(1-m[:α])),
                                  description = "GDP", tex_label = "y_*")
        m <= SteadyStateParameter(:xstar, (1-(1-m[:δ])*exp(-m[:γ]))*m[:kstar],
                                  description = "Investment", tex_label = "x_*")
        m <= SteadyStateParameter(:cstar, (m.grids[:weights_total].*het[:μstar].value)'*(min.(1.0 ./ het[:lstar].value, repeat(m.grids[:xgrid].points, ns) .+ m[:η])),
                                  description = "Steady-state consumption", tex_label = "c_*")
    end
    m <= SteadyStateParameter(:lstar, 1/m[:cstar],
                              description = "Steady-state expected discounted
                              marginal utility of consumption", tex_label = "l_*")
    m <= SteadyStateParameter(:bg, m[:BoverY]*m[:ystar]*exp(m[:γ]), description = "Govt Debt",
                              tex_label = "bg")
    m <= SteadyStateParameter(:Tg, (exp(-m[:γ]) - 1/(1+m[:r]))*m[:bg] + (1-(1/m[:g]))*m[:ystar],
                              description = "Net lump sump taxes", tex_label = "Tg")
    m <= SteadyStateParameter(:Tstar, m[:Rkstar]*m[:kstar]*exp(-m[:γ]) - m[:xstar] - m[:Tg],
                              description = "Net transfer to households", tex_label = "T_*")

    return m
end
