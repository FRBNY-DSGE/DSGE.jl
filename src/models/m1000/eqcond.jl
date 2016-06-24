"""
```
eqcond(m::Model990)
```

Expresses the equilibrium conditions in canonical form using Γ0, Γ1, C, Ψ, and Π matrices.
Using the mappings of states/equations to integers defined in m990.jl, coefficients are
specified in their proper positions.

### Outputs

* `Γ0` (`n_states` x `n_states`) holds coefficients of current time states.
* `Γ1` (`n_states` x `n_states`) holds coefficients of lagged states.
* `C`  (`n_states` x `1`) is a vector of constants
* `Ψ`  (`n_states` x `n_shocks_exogenous`) holds coefficients of iid shocks.
* `Π`  (`n_states` x `n_states_expectational`) holds coefficients of expectational states.
"""
function eqcond(m::Model990)
    endo = m.endogenous_states
    exo  = m.exogenous_shocks
    ex   = m.expected_shocks
    eq   = m.equilibrium_conditions
    
    #SUMMARY
    neq = 8
    neta = 2
    neps = 3

    Γ0 = zeros(neq,neq)
    Γ1 = zeros(neq,neq)
    C  = zeros(neq, 1)
    Ψ  = zeros(neq,neps)
    Π  = zeros(neq,neta)

    ### ENDOGENOUS STATES ###

    ### 1. Consumption Euler Equation

    Γ0[eq[:eq_euler], endo[:y_t]]  = 1.
    Γ0[eq[:eq_euler], endo[:R_t]] = 1/tau
    Γ0[eq[:eq_euler], endo[:g_t]] = -(1-rho_g)
    Γ0[eq[:eq_euler], endo[:z_t]] = -rho_z/tau
    Γ0[eq[:eq_euler], endo[:Ey_t1]] = -1
    Γ0[eq[:eq_euler], endo[:Epi_t1]] = -1/tau

    ### 2. NK Phillips Curve

    # Sticky prices and wages
    Γ0[eq[:nk_pcurve], endo[:y_t]] = -kappa
    Γ0[eq[:nk_pcurve], endo[:pi_t]] = 1
    Γ0[eq[:nk_pcurve], endo[:g_t]] = kappa
    Γ0[eq[:nk_pcurve], endo[:Epi_t1]] = -bet

    ### 3. Monetary Policy Rule

    Γ0[eq[:mp_rule], endo[:y_t]] = -(1-rho_R)*psi2
    Γ0[eq[:mp_rule], endo[:pi_t]] = -(1-rho_R)*psi1
    Γ0[eq[:mp_rule], endo[:R_t]] = 1
    Γ0[eq[:mp_rule], endo[:g_t]] = (1-rho_R)*psi2
    Γ1[eq[:mp_rule], endo[:R_t]] = rho_R
    Ψ[eq[:mp_rule], exo[:R_sh]] = 1

    ### 4. Shock 1

    Γ0[eq[:shock_1], endo[:y1_t]] = 1
    Γ1[eq[:shock_1], endo[:y_t]] = 1

    ### 5. Shock 2

    Γ0[eq[:shock_2], endo[:g_t]] = 1
    Γ1[eq[:shock_2], endo[:g_t]] = rho_g
    Ψ[eq[:shock_2], exo[:g_sh]] = 1

    ### 6. Shock 3

    Γ0[eq[:shock_3], endo[:z_t]] = 1
    Γ1[eq[:shock_3], endo[:z_t]] = rho_z
    Ψ[eq[:shock_3], exo[:z_sh]] = 1

    ### 7. Shock 4

    Γ0[eq[:shock_4], endo[:y_t]] = 1
    Γ1[eq[:shock_4], endo[:Ey_t1]] = 1
    Π[eq[:shock_4], ex[:ey_sh]] = 1

    ### 8. Shock 5

    Γ0[eq[:shock_5], endo[:pi_t]] = 1
    Γ1[eq[:shock_5], endo[:Epi_t1]] = 1
    Π[eq[:shock_5], ex[:epi_sh]] = 1

    return Γ0, Γ1, C, Ψ, Π

end
