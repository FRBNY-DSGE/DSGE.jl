"""
```
eqcond(m::Schorf)
```

Expresses the equilibrium conditions in canonical form using Γ0, Γ1, C, Ψ, and Π matrices.
Using the mappings of states/equations to integers defined in schorf.jl, coefficients are
specified in their proper positions.

### Outputs

* `Γ0` (`n_states` x `n_states`) holds coefficients of current time states.
* `Γ1` (`n_states` x `n_states`) holds coefficients of lagged states.
* `C`  (`n_states` x `1`) is a vector of constants
* `Ψ`  (`n_states` x `n_shocks_exogenous`) holds coefficients of iid shocks.
* `Π`  (`n_states` x `n_states_expectational`) holds coefficients of expectational states.
"""
function eqcond(m::AnSchorfheide)
    endo = m.endogenous_states
    exo  = m.exogenous_shocks
    ex   = m.expected_shocks
    eq   = m.equilibrium_conditions

    Γ0 = zeros(neq,neq)
    Γ1 = zeros(neq,neq)
    C  = zeros(neq, 1)
    Ψ  = zeros(neq,neps)
    Π  = zeros(neq,neta)

    ### ENDOGENOUS STATES ###

    ### 1. Consumption Euler Equation

    Γ0[eq[:eq_euler], endo[:y_t]]  = 1.
    Γ0[eq[:eq_euler], endo[:R_t]] = 1/m[:τ]
    Γ0[eq[:eq_euler], endo[:g_t]] = -(1-m[:ρ_g])
    Γ0[eq[:eq_euler], endo[:z_t]] = -m[:ρ_z]/m[:τ]
    Γ0[eq[:eq_euler], endo[:Ey_t1]] = -1
    Γ0[eq[:eq_euler], endo[:Eπ_t1]] = -1/m[:τ]

    ### 2. NK Phillips Curve

    Γ0[eq[:eq_phillips], endo[:y_t]] = -m[:κ]
    Γ0[eq[:eq_phillips], endo[:π_t]] = 1
    Γ0[eq[:eq_phillips], endo[:g_t]] = m[:κ]
    Γ0[eq[:eq_phillips], endo[:Eπ_t1]] = -1/(1+m[:rA]/400)

    ### 3. Monetary Policy Rule

    Γ0[eq[:eq_mp], endo[:y_t]] = -(1-m[:ρ_R])*m[:ψ_2]
    Γ0[eq[:eq_mp], endo[:π_t]] = -(1-m[:ρ_R])*m[:ψ_1]
    Γ0[eq[:eq_mp], endo[:R_t]] = 1
    Γ0[eq[:eq_mp], endo[:g_t]] = (1-m[:ρ_R])*m[:ψ_2]
    Γ1[eq[:eq_mp], endo[:R_t]] = m[:ρ_R]
    Ψ[eq[:eq_mp], exo[:R_sh]] = 1

    ### 4. Output lag

    Γ0[eq[:eq_y_t1], endo[:y_t1]] = 1
    Γ1[eq[:eq_y_t1], endo[:y_t]] = 1

    ### 5. Government spending

    Γ0[eq[:eq_g], endo[:g_t]] = 1
    Γ1[eq[:eq_g], endo[:g_t]] = m[:ρ_g]
    Ψ[eq[:eq_g], exo[:g_sh]] = 1

    ### 6. Technology

    Γ0[eq[:eq_z], endo[:z_t]] = 1
    Γ1[eq[:eq_z], endo[:z_t]] = m[:ρ_z]
    Ψ[eq[:eq_z], exo[:z_sh]] = 1

    ### 7. Expected output

    Γ0[eq[:eq_Ey], endo[:y_t]] = 1
    Γ1[eq[:eq_Ey], endo[:Ey_t1]] = 1
    Π[eq[:eq_Ey], ex[:Ey_sh]] = 1

    ### 8. Expected inflation

    Γ0[eq[:eq_Eπ], endo[:π_t]] = 1
    Γ1[eq[:eq_Eπ], endo[:Eπ_t1]] = 1
    Π[eq[:eq_Eπ], ex[:Eπ_sh]] = 1

    return Γ0, Γ1, C, Ψ, Π
end
