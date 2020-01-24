"""
```
pseudo_measurement(m::SmetsWoutersOrig{T},
    TTT::Matrix{T}, RRR::Matrix{T}, CCC::Vector{T}) where {T<:AbstractFloat}
```

Assign pseudo-measurement equation (a linear combination of states):

```
x_t = ZZ_pseudo*s_t + DD_pseudo
```
"""
function pseudo_measurement(m::SmetsWoutersOrig{T},
                            TTT::Matrix{T},
                            RRR::Matrix{T},
                            CCC::Vector{T}) where {T<:AbstractFloat}

    endo      = m.endogenous_states
    endo_addl = m.endogenous_states_augmented
    pseudo    = m.pseudo_observables

    _n_states = n_states_augmented(m)
    _n_pseudo = n_pseudo_observables(m)

    # Initialize pseudo ZZ and DD matrices
    ZZ_pseudo = zeros(_n_pseudo, _n_states)
    DD_pseudo = zeros(_n_pseudo)

    ## Output
    ZZ_pseudo[pseudo[:y_t],endo[:y_t]] = 1.

    ## Marginal Cost
    ZZ_pseudo[pseudo[:MarginalCost],endo[:mc_t]] = 1.

    ## Wages
    ZZ_pseudo[pseudo[:Wages],endo[:w_t]] = 1.

    ## Hours
    ZZ_pseudo[pseudo[:Hours],endo[:L_t]] = 1.

    ## z_t
    ZZ_pseudo[pseudo[:z_t], endo[:z_t]] = 1.

    ## Nominal FFR
    ZZ_pseudo[pseudo[:NominalFFR], endo[:R_t]] = 1.
    DD_pseudo[pseudo[:NominalFFR]] = m[:Rstarn]

    ## Nominal Wage Growth
    ZZ_pseudo[pseudo[:NominalWageGrowth],endo[:w_t]] = 1.
    ZZ_pseudo[pseudo[:NominalWageGrowth],endo_addl[:w_t1]] = -1.
    ZZ_pseudo[pseudo[:NominalWageGrowth],endo[:π_t]] = 1.
    DD_pseudo[pseudo[:NominalWageGrowth]]            = m[:π_star].value + m[:γ].value

    ZZ_pseudo[pseudo[:LaborProductivityGrowthNoME], endo[:y_t]]       = 1.
    ZZ_pseudo[pseudo[:LaborProductivityGrowthNoME], endo_addl[:y_t1]] = -1.
    ZZ_pseudo[pseudo[:LaborProductivityGrowthNoME], endo[:L_t]]       = -1
    ZZ_pseudo[pseudo[:LaborProductivityGrowthNoME], endo_addl[:L_t1]] = 1.
    DD_pseudo[pseudo[:LaborProductivityGrowthNoME]]                   = m[:γ].value

    ZZ_pseudo[pseudo[:laborshare_t], endo[:w_t]] = 1.
    ZZ_pseudo[pseudo[:laborshare_t], endo[:L_t]] = 1.
    ZZ_pseudo[pseudo[:laborshare_t], endo[:y_t]] = -1.
    # wl_c = wstar^h * Lstar / cstar = 1 / λ_w * (wstar * Lstar / cstar)
    # c_y  = cstar / ystar
    # wstar * Lstar / cstar = wl_c * λ_w
    # wstar * Lstar / cstar * (cstar / ystar) = wstar * Lstar / ystar
    DD_pseudo[pseudo[:laborshare_t]] = 100. * log(m[:λ_w] * m[:wl_c] * m[:c_y])

    return PseudoMeasurement(ZZ_pseudo, DD_pseudo)
end

function pseudo_measurement(m::SmetsWoutersOrig{T},
                            TTTs::Vector{Matrix{T}},
                            RRRs::Vector{Matrix{T}},
                            CCCs::Vector{Vector{T}}) where {T<:AbstractFloat}
    return [pseudo_measurement(m, TTTs[1], RRRs[1], CCCs[1]),
                               pseudo_measurement(m, TTTs[2], RRRs[2], CCCs[2])]
end
