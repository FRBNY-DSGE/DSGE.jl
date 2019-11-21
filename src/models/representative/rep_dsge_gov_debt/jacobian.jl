function jacobian(m::RepDSGEGovDebt)
    endo = augment_model_states(m.endogenous_states, n_model_states(m))
    eq   = m.equilibrium_conditions

    # Load in parameters/steady-state parameters
    c_scalar   = m[:cstar].value
    ell_scalar = m[:lstar].value
    ω          = m[:ωstar].value
    H          = m[:H].value
    ϕh         = m[:ϕh].value
    lam        = 1/c_scalar
    ϕ          = lam*ω/(H^ϕh)
    βrank      = exp(m[:γ].value)/(1+m[:r])
    y          = m[:ystar].value
    g          = m[:g].value
    Rk         = m[:Rkstar].value
    T          = m[:Tstar].value
    k          = m[:kstar].value
    x          = m[:xstar].value
    Tg         = m[:Tg].value
    spp        = m[:spp].value
    R          = m[:r].scaledvalue + 1
    γ          = m[:γ].scaledvalue
    δ          = m[:δ].value
    κw         = m[:κ_w].value
    κp         = m[:κ_p].value
    α          = m[:α].value
    ρR         = m[:ρR].value
    ψπ         = m[:ψπ].value
    ψy         = m[:ψy].value
    δb         = m[:δb].value
    bg         = m[:bg].value
    ρB         = m[:ρ_b].value
    ρG         = m[:ρ_g].value
    ρZ         = m[:ρ_z].value
    ρμ         = m[:ρ_μ].value
    ρlamw      = m[:ρ_λ_w].value
    ρlamf      = m[:ρ_λ_f].value
    ρmon       = m[:ρ_rm].value

    # make the Jacobian
    JJ = zeros(n_model_states(m), 2*n_model_states(m))

    # Euler equation
    JJ[eq[:eq_euler],endo[:l′_t]] = 1.
    JJ[eq[:eq_euler],endo[:z′_t]] = -1.
    JJ[eq[:eq_euler],endo[:b_t]]  = ell_scalar
    JJ[eq[:eq_euler],endo[:l_t]]  = -1.
    JJ[eq[:eq_euler],endo[:R_t]]  = 1.

    # mkt clearing
    JJ[eq[:eq_market_clearing],endo[:y_t]]   = y/g
    JJ[eq[:eq_market_clearing],endo[:g_t]]   = -y/g
    JJ[eq[:eq_market_clearing],endo[:I_t]]   = -x
    JJ[eq[:eq_market_clearing],endo[:l_t]] = c_scalar

    # lambda = average marginal utility
    JJ[eq[:eq_lambda],endo[:margutil_t]] = lam
    JJ[eq[:eq_lambda],endo[:l_t]] = -1/c_scalar

    # transfer
    JJ[eq[:eq_transfers],endo[:t_t]]  = T
    JJ[eq[:eq_transfers],endo[:capreturn_t]]  = -Rk*k
    JJ[eq[:eq_transfers],endo[:k_t]]  = -Rk*k
    JJ[eq[:eq_transfers],endo[:z_t]]   = Rk*k
    JJ[eq[:eq_transfers],endo[:I_t]]   = x
    JJ[eq[:eq_transfers],endo[:mc_t]]  = y
    JJ[eq[:eq_transfers],endo[:tg_t]]  = Tg

    # investment
    JJ[eq[:eq_investment],endo[:Q_t]]  = 1.
    JJ[eq[:eq_investment],endo[:μ_t]]  = 1.
    JJ[eq[:eq_investment],endo[:I′_t]] = spp*(exp(3*γ))/R
    JJ[eq[:eq_investment],endo[:z′_t]] = spp*(exp(3*γ))/R
    JJ[eq[:eq_investment],endo[:I_t]]  = -spp*(exp(3*γ))/R - spp*exp(2*γ)
    JJ[eq[:eq_investment],endo[:I_t1]] = spp*exp(2*γ)
    JJ[eq[:eq_investment],endo[:z_t]]  = -spp*exp(2*γ)

    # tobin's q
    JJ[eq[:eq_tobin_q],endo[:margutil_t]]  = 1.
    JJ[eq[:eq_tobin_q],endo[:margutil′_t]] = -1.
    JJ[eq[:eq_tobin_q],endo[:Q_t]]    = 1.
    JJ[eq[:eq_tobin_q],endo[:z′_t]]   = 1.
    JJ[eq[:eq_tobin_q],endo[:capreturn′_t]]  = -Rk/R
    JJ[eq[:eq_tobin_q],endo[:Q′_t]]   = -(1-δ)/R

    # capital accumulation
    JJ[eq[:eq_capital_accumulation],endo[:k′_t]] = 1.
    JJ[eq[:eq_capital_accumulation],endo[:k_t]]  = -(1-δ)
    JJ[eq[:eq_capital_accumulation],endo[:z_t]]   = (1-δ)
    JJ[eq[:eq_capital_accumulation],endo[:μ_t]]  = -x/k
    JJ[eq[:eq_capital_accumulation],endo[:I_t]]   = -x/k

    # wage phillips curve
    JJ[eq[:eq_wage_phillips],endo[:π_w_t]]  = -1.
    JJ[eq[:eq_wage_phillips],endo[:λ_w_t]] = 1.
    JJ[eq[:eq_wage_phillips],endo[:L_t]]   = κw*ϕh
    JJ[eq[:eq_wage_phillips],endo[:w_t]]    = -κw
    JJ[eq[:eq_wage_phillips],endo[:π_w′_t]]  = βrank
    JJ[eq[:eq_wage_phillips],endo[:margutil_t]]  = -κw

    # price phillips curve
    JJ[eq[:eq_price_phillips],endo[:π_t]]   = -1.
    JJ[eq[:eq_price_phillips],endo[:mc_t]]   = κp
    JJ[eq[:eq_price_phillips],endo[:λ_f_t]] = 1.
    JJ[eq[:eq_price_phillips],endo[:π′_t]]  = 1/R

    # marginal cost
    JJ[eq[:eq_marginal_cost],endo[:mc_t]] = 1.
    JJ[eq[:eq_marginal_cost],endo[:w_t]]  = -(1-α)
    JJ[eq[:eq_marginal_cost],endo[:capreturn_t]] = -α

    # gdp
    JJ[eq[:eq_gdp],endo[:y_t]]  = 1.
    JJ[eq[:eq_gdp],endo[:z_t]]  = α
    JJ[eq[:eq_gdp],endo[:k_t]] = -α
    JJ[eq[:eq_gdp],endo[:L_t]] = -(1-α)

    # optimal k/l ratio
    JJ[eq[:eq_optimal_kl],endo[:capreturn_t]] = 1.
    JJ[eq[:eq_optimal_kl],endo[:w_t]]  = -1.
    JJ[eq[:eq_optimal_kl],endo[:L_t]] = -1.
    JJ[eq[:eq_optimal_kl],endo[:k_t]] = 1.
    JJ[eq[:eq_optimal_kl],endo[:z_t]]  = -1.

    # taylor rule
    JJ[eq[:eq_taylor],endo[:i_t]]   = -1.
    JJ[eq[:eq_taylor],endo[:i_t1]]  = ρR
    JJ[eq[:eq_taylor],endo[:π_t]]   = (1-ρR)*ψπ
    JJ[eq[:eq_taylor],endo[:y_t]]    = (1-ρR)*ψy
    JJ[eq[:eq_taylor],endo[:y_t1]]  = -(1-ρR)*ψy
    JJ[eq[:eq_taylor],endo[:z_t]]    = (1-ρR)*ψy
    JJ[eq[:eq_taylor],endo[:rm_t]] = 1.

    # fisher eqn
    JJ[eq[:eq_fisher],endo[:R_t]]  = 1.
    JJ[eq[:eq_fisher],endo[:π′_t]] = 1.
    JJ[eq[:eq_fisher],endo[:i_t]]  = -1.

    # wage inflation
    JJ[eq[:eq_nominal_wage_inflation],endo[:π_w_t]] = 1.
    JJ[eq[:eq_nominal_wage_inflation],endo[:π_t]]  = -1.
    JJ[eq[:eq_nominal_wage_inflation],endo[:z_t]]   = -1.
    JJ[eq[:eq_nominal_wage_inflation],endo[:w_t]]   = -1.
    JJ[eq[:eq_nominal_wage_inflation],endo[:w_t1]]  = 1.

    # fiscal rule
    JJ[eq[:eq_fiscal_rule],endo[:tg_t]]  = -Tg
    JJ[eq[:eq_fiscal_rule],endo[:R_t]]   = δb*bg/R
    JJ[eq[:eq_fiscal_rule],endo[:bg_t]]  = δb*bg*ℯ^(-γ)
    JJ[eq[:eq_fiscal_rule],endo[:z_t]]   = -δb*bg*ℯ^(-γ)
    JJ[eq[:eq_fiscal_rule],endo[:y_t]]   = δb*(1-(1/g))*y
    JJ[eq[:eq_fiscal_rule],endo[:g_t]]   = δb*y/g

    # govt budget constraint
    JJ[eq[:eq_g_budget_constraint],endo[:bg′_t]] = -(bg/R)
    JJ[eq[:eq_g_budget_constraint],endo[:R_t]]   = bg/R
    JJ[eq[:eq_g_budget_constraint],endo[:bg_t]]  = bg*ℯ^(-γ)
    JJ[eq[:eq_g_budget_constraint],endo[:z_t]]   = -bg*ℯ^(-γ)
    JJ[eq[:eq_g_budget_constraint],endo[:y_t]]   = (1-(1/g))*y
    JJ[eq[:eq_g_budget_constraint],endo[:g_t]]   = y/g
    JJ[eq[:eq_g_budget_constraint],endo[:tg_t]]  = -Tg

    # update lagged variables
    JJ[eq[:LR],endo[:R′_t1]] = 1.
    JJ[eq[:LR],endo[:R_t]]   = -1.

    JJ[eq[:LI],endo[:i′_t1]] = 1.
    JJ[eq[:LI],endo[:i_t]]   = -1.

    JJ[eq[:LY],endo[:y′_t1]] = 1.
    JJ[eq[:LY],endo[:y_t]]   = -1.

    JJ[eq[:LW],endo[:w′_t1]] = 1.
    JJ[eq[:LW],endo[:w_t]]   = -1.

    JJ[eq[:LX],endo[:I′_t1]] = 1.
    JJ[eq[:LX],endo[:I_t]]   = -1.

    # discount factor shock
    JJ[eq[:eq_b],endo[:b′_t]] = 1.
    JJ[eq[:eq_b],endo[:b_t]]  = -ρB

    # g/y shock
    JJ[eq[:eq_g],endo[:g′_t]] = 1.
    JJ[eq[:eq_g],endo[:g_t]]  = -ρG

    # tfp growth shock
    JJ[eq[:eq_z],endo[:z′_t]] = 1.
    JJ[eq[:eq_z],endo[:z_t]]  = -ρZ

    # investment shock
    JJ[eq[:eq_μ],endo[:μ′_t]] = 1.
    JJ[eq[:eq_μ],endo[:μ_t]]  = -ρμ

    # wage mkup shock
    JJ[eq[:eq_λ_w],endo[:λ_w′_t]] = 1.
    JJ[eq[:eq_λ_w],endo[:λ_w_t]]  = -ρlamw

    # price mkup shock
    JJ[eq[:eq_λ_f],endo[:λ_f′_t]] = 1.
    JJ[eq[:eq_λ_f],endo[:λ_f_t]] = -ρlamf

    # monetary policy shock
    JJ[eq[:eq_rm],endo[:rm′_t]] = 1.
    JJ[eq[:eq_rm],endo[:rm_t]]  = -ρmon

    if !m.testing
        JJ  = normalize(m, JJ)
    end

    return JJ
end

function compose_normalization_matrices(m::RepDSGEGovDebt)
    Qleft     = eye(n_model_states(m))
    Qx        = eye(n_backward_looking_states(m))
    Qy        = eye(n_jumps(m))

    Qright    = cat(Qx', Qy', Qx', Qy', dims = [1,2])

    return Qx, Qy, Qleft, Qright
end

function normalize(m::RepDSGEGovDebt, JJ::Matrix{Float64})
    Qx, _, Qleft, Qright = compose_normalization_matrices(m)

    m <= Setting(:n_predetermined_variables, size(Qx, 1))

    Jac1 = Qleft*JJ*Qright

    return Jac1
end
