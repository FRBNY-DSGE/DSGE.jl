"""
```
eqcond(m::Schorf)
```

Expresses the equilibrium conditions in canonical form using Γ0, Γ1, C, Ψ, and Π matrices.
Using the mappings of states/equations to integers defined in an_schorfheide.jl, coefficients are
specified in their proper positions.

### Outputs

* `Γ0` (`n_states` x `n_states`) holds coefficients of current time states.
* `Γ1` (`n_states` x `n_states`) holds coefficients of lagged states.
* `C`  (`n_states` x `1`) is a vector of constants
* `Ψ`  (`n_states` x `n_shocks_exogenous`) holds coefficients of iid shocks.
* `Π`  (`n_states` x `n_states_expectational`) holds coefficients of expectational states.
"""
function eqcond(m::AnSchorfheide) # do not edit this function definition
    return eqcond(m, 1)
end

function eqcond(m::AnSchorfheide, reg::Int) # do not edit these inputs
    endo = m.endogenous_states # for convenience, make references to these fields
    exo  = m.exogenous_shocks  # which hold indices, e.g. endo[:y_t] will return an Int
    ex   = m.expected_shocks   # for the index of the `y_t` variable.
    eq   = m.equilibrium_conditions # you can rename these abbreviations if you want

    Γ0 = zeros(n_states(m), n_states(m)) # initialize containers (do not edit)
    Γ1 = zeros(n_states(m), n_states(m))
    C  = zeros(n_states(m))
    Ψ  = zeros(n_states(m), n_shocks_exogenous(m))
    Π  = zeros(n_states(m), n_shocks_expectational(m))

    # Regime-switching parameters (you should probably not edit this call)
    for para in m.parameters      # if you're new to DSGE.jl, then revisit this block of code later.
        if !isempty(para.regimes) # if you don't need regime-switching, then you can comment out this block of code
            if (haskey(get_settings(m), :model2para_regime) ? haskey(get_setting(m, :model2para_regime), para.key) : false)
                toggle_regime!(para, reg, get_setting(m, :model2para_regime)[para.key])
            else
                toggle_regime!(para, reg)
            end
        end
    end

    # You should definitely make edits to the following code!

    ### ENDOGENOUS STATES ###

    ### 1. Consumption Euler Equation

    Γ0[eq[:eq_euler], endo[:y_t]]  = 1.
    Γ0[eq[:eq_euler], endo[:R_t]] = 1/m[:τ]
    Γ0[eq[:eq_euler], endo[:g_t]] = -(1-m[:ρ_g])
    Γ0[eq[:eq_euler], endo[:z_t]] = -m[:ρ_z]/m[:τ]
    Γ0[eq[:eq_euler], endo[:Ey_t]] = -1
    Γ0[eq[:eq_euler], endo[:Eπ_t]] = -1/m[:τ]

    ### 2. NK Phillips Curve

    Γ0[eq[:eq_phillips], endo[:y_t]] = -m[:κ]
    Γ0[eq[:eq_phillips], endo[:π_t]] = 1
    Γ0[eq[:eq_phillips], endo[:g_t]] = m[:κ]
    Γ0[eq[:eq_phillips], endo[:Eπ_t]] = -1/(1+m[:rA]/400)

    ### 3. Monetary Policy Rule

    Γ0[eq[:eq_mp], endo[:y_t]] = -(1-m[:ρ_R])*m[:ψ_2]
    Γ0[eq[:eq_mp], endo[:π_t]] = -(1-m[:ρ_R])*m[:ψ_1]
    Γ0[eq[:eq_mp], endo[:R_t]] = 1
    Γ0[eq[:eq_mp], endo[:g_t]] = (1-m[:ρ_R])*m[:ψ_2]
    Γ1[eq[:eq_mp], endo[:R_t]] = m[:ρ_R]
    Ψ[eq[:eq_mp], exo[:rm_sh]] = 1

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
    Γ1[eq[:eq_Ey], endo[:Ey_t]] = 1
    Π[eq[:eq_Ey], ex[:Ey_sh]] = 1

    ### 8. Expected inflation

    Γ0[eq[:eq_Eπ], endo[:π_t]] = 1
    Γ1[eq[:eq_Eπ], endo[:Eπ_t]] = 1
    Π[eq[:eq_Eπ], ex[:Eπ_sh]] = 1

    # Ensure parameter regimes are in 1 at the end (you should probably not edit)
    for para in m.parameters      # if you're new to DSGE.jl, then revisit this block of code later.
        if !isempty(para.regimes) # if you don't need regime-switching, then you can comment out this block of code
            toggle_regime!(para, 1)
        end
    end

    return Γ0, Γ1, C, Ψ, Π # do not edit this return
end
