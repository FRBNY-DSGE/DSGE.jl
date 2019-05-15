# Lightweight Representative Agent Version of HetDSGEGovDebt
mutable struct RepDSGEGovDebt{T} <: AbstractRepModel{T}
    parameters::ParameterVector{T}               # vector of all time-invariant model parameters
    steady_state::ParameterVector{T}             # model steady-state values

    grids::OrderedDict{Symbol,Union{Grid, Array}}
    keys::OrderedDict{Symbol,Int}                    # Human-readable names for all the model
                                                     # parameters and steady-states

    state_variables::Vector{Symbol}                  # Vector of symbols of the state variables
    jump_variables::Vector{Symbol}                   # Vector of symbols of the jump variables

    # Vector of unnormalized ranges of indices
    endogenous_states::OrderedDict{Symbol,Int}
    exogenous_shocks::OrderedDict{Symbol,Int}
    expected_shocks::OrderedDict{Symbol,Int}
    equilibrium_conditions::OrderedDict{Symbol,Int}
    endogenous_states_augmented::OrderedDict{Symbol, Int}
    observables::OrderedDict{Symbol,Int}

    spec::String                                     # Model specification number (eg "m990")
    subspec::String                                  # Model subspecification (eg "ss0")
    settings::Dict{Symbol,Setting}                   # Settings/flags for computation
    testing::Bool                                    # Whether we are in testing mode or not
    observable_mappings::OrderedDict{Symbol, Observable}
end

description(m::RepDSGEGovDebt) = "RepDSGEGovDebt, $(m.subspec)"

function RepDSGEGovDebt(het::HetDSGEGovDebt)

    @assert !isnan(het[:βstar].value) "Must have already called steadystate! on `het` before instantiating a RepDSGEGovDebt model"

    # Model-specific specifications
    spec               = "rep_dsge"
    subspec            = het.subspec
    settings           = het.settings

    # initialize empty model
    m = RepDSGEGovDebt{Float64}(
            # model parameters and steady state values
            Vector{AbstractParameter{Float64}}(), Vector{Float64}(),
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
            het.testing,
            het.observable_mappings)

    # Set settings
    default_settings!(m)
    model_settings!(m, het)

    # Initialize model indices
    init_model_indices!(m)

    # Initialize parameters
    init_parameters!(m, het)

    # Initialize grids
    init_grids!(m, het)

    # # Solve for the steady state
    # steadystate!(m, het)

    return m
end

function model_settings!(m::RepDSGEGovDebt, het::HetDSGEGovDebt)
    m.settings = het.settings

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
end

function init_model_indices!(m::RepDSGEGovDebt)

    # Endogenous states
    states = collect([:k′_t, :R′_t1, :i′_t1, :y′_t1, :w′_t1, :I′_t1, :bg′_t,
                      :b′_t, :g′_t, :z′_t, :μ′_t, :λ_w′_t, :λ_f′_t, :rm′_t])

    jumps = collect([:R′_t, :i′_t, :t′_t, :w′_t, :L′_t, :π′_t, :π_w′_t, :mu′_t, :y′_t, :I′_t,
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
    endogenous_states_augmented = [:i_t1, :c_t, :c_t1]

    # Observables
    observables = keys(m.observable_mappings)

    ########################################################################################
    # Setting indices of endogenous_states and equilibrium conditions manually for now

    # ATTN: Probably unnecessary. Can do it the enumerate way (below) since it's RANK
    # Delete when comfortable
    # setup_indices!(m)
    endo   = m.endogenous_states
    eqcond = equilibrium_conditions
    ########################################################################################

    for (i,k) in enumerate(endogenous_states);m.endogenous_states[k] = i end
    for (i,k) in enumerate(endogenous_states_augmented); m.endogenous_states_augmented[k] = i+length(endogenous_states) end
    for (i,k) in enumerate(exogenous_shocks); m.exogenous_shocks[k]  = i end
    for (i,k) in enumerate(observables);      m.observables[k]       = i end
    for (i,k) in enumerate(equilibrium_conditions);      m.equilibrium_conditions[k]      = i end
end

# ATTN: Probably unneeded
function setup_indices!(m::RepDSGEGovDebt)
    nx = get_setting(m, :nx)
    ns = get_setting(m, :ns)
    endo = m.endogenous_states
    eqconds = m.equilibrium_conditions

    # endogenous scalar-valued states
    endo[:k′_t]   = 1             # capital –dont get confused with steadystate object K
    endo[:R′_t1]  = 2             # lagged real interest rate
    endo[:i′_t1]  = 3             # lagged nominal interest rate
    endo[:y′_t1]  = 4             # lagged gdp
    endo[:w′_t1]  = 5             # lag real wages
    endo[:I′_t1]  = 6             # lag investment–don't get this confused with i, the nominal interest rate
    endo[:bg′_t]  = 7             # govt debt

    # exogenous scalar-valued states:
    endo[:b′_t]   = 8             # discount factor shock
    endo[:g′_t]   = 9             # govt spending
    endo[:z′_t]   = 10            # tfp growth
    endo[:μ′_t]   = 11            # investment shock
    endo[:λ_w′_t] = 12            # wage markup
    endo[:λ_f′_t] = 13            # price markup
    endo[:rm′_t]  = 14            # monetary policy shock
    # endo[:c′_t1]  = 14            # lagged consumption

    # scalar-valued jumps
    endo[:R′_t]   = 15            # real interest rate
    endo[:i′_t]   = 16            # nominal interest rate
    endo[:t′_t]   = 17            # transfers + dividends
    endo[:w′_t]   = 18            # real wage
    endo[:L′_t]   = 19            # hours worked
    endo[:π′_t]   = 20            # inflation
    endo[:π_w′_t] = 21            # nominal wage inflation
    endo[:mu′_t]  = 22            # average marginal utility
    endo[:y′_t]   = 23            # gdp
    endo[:I′_t]   = 24            # investment
    endo[:mc′_t]  = 25            # marginal cost - this is ζ in HetDSGEGovDebtₖd.pdf
    endo[:Q′_t]   = 26            # Tobin's qfunction
    endo[:capreturn′_t] = 27      # return on capital
    endo[:l′_t]   = 28            # marginal utility, now a scalar
    endo[:tg′_t]  = 29 # lump sum tax

    eqconds[:eq_euler]                = 1
    eqconds[:eq_market_clearing]      = 2
    eqconds[:eq_lambda]               = 3
    eqconds[:eq_transfers]            = 4
    eqconds[:eq_investment]           = 5
    eqconds[:eq_tobin_q]              = 6
    eqconds[:eq_capital_accumulation] = 7
    eqconds[:eq_wage_phillips]        = 8
    eqconds[:eq_price_phillips]       = 9
    eqconds[:eq_marginal_cost]        = 10
    eqconds[:eq_gdp]                  = 11
    eqconds[:eq_optimal_kl]           = 12
    eqconds[:eq_taylor]               = 13
    eqconds[:eq_fisher]               = 14
    eqconds[:eq_nominal_wage_inflation] = 15
    eqconds[:eq_fiscal_rule]          = 16
    eqconds[:eq_g_budget_constraint]  = 17

    # lagged variables
    eqconds[:LR] = 18
    eqconds[:LI] = 19
    eqconds[:LY] = 20
    eqconds[:LW] = 21
    eqconds[:LX] = 22
    # shocks
    eqconds[:eq_b]   = 23
    eqconds[:eq_g]   = 24
    eqconds[:eq_z]   = 25
    eqconds[:eq_μ]   = 26
    eqconds[:eq_λ_w] = 27
    eqconds[:eq_λ_f] = 28
    eqconds[:eq_rm]  = 29

    m.endogenous_states = deepcopy(endo)
end

function init_parameters!(m::RepDSGEGovDebt, het::HetDSGEGovDebt)
    for param in het.parameters
        m <= param
    end
end

function init_grids!(m::RepDSGEGovDebt, het::HetDSGEGovDebt)
    m.grids = het.grids
end

function steadystate!(m::RepDSGEGovDebt, het::HetDSGEGovDebt)
    ns = get_setting(m, :ns)

    m <= SteadyStateParameter(:cstar, (m.grids[:weights_total].*het[:μstar].value)'*(min.(1./het[:lstar].value, repeat(m.grids[:xgrid].points, ns) .+ m[:η])),
                              description = "Steady-state consumption", tex_label = "c_*")
    m <= SteadyStateParameter(:lstar, 1/m[:cstar],
                              description = "Steady-state expected discounted
                              marginal utility of consumption", tex_label = "l_*")
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
    m <= SteadyStateParameter(:xstar, (1-(1-m[:δ])*exp(-m[:γ]))*m[:kstar],
                              description = "Investment", tex_label = "x_*")
    m <= SteadyStateParameter(:ystar, (exp(-m[:α]*m[:γ])*(m[:kstar]^m[:α])*m[:H]^(1-m[:α])),
                              description = "GDP", tex_label = "y_*")
    m <= SteadyStateParameter(:bg, m[:BoverY]*m[:ystar]*exp(m[:γ]), description = "Govt Debt",
                              tex_label = "bg")
    m <= SteadyStateParameter(:Tg, (exp(-m[:γ]) - 1/(1+m[:r]))*m[:bg] + (1-(1/m[:g]))*m[:ystar],
                              description = "Net lump sump taxes", tex_label = "Tg")
    m <= SteadyStateParameter(:Tstar, m[:Rkstar]*m[:kstar]*exp(-m[:γ]) - m[:xstar] - m[:Tg],
                              description = "Net transfer to households", tex_label = "T_*")
    return m
end
