"""
```
pseudo_measurement(m::AnSchorfheide{T},
    TTT::Matrix{T}, RRR::Matrix{T}, CCC::Vector{T}) where {T<:AbstractFloat}
```

Assign pseudo-measurement equation (a linear combination of states):

```
x_t = ZZ_pseudo*s_t + DD_pseudo
```
"""
function pseudo_measurement(m::AnSchorfheide{T},
                            TTT::Matrix{T},
                            RRR::Matrix{T},
                            CCC::Vector{T}; reg::Int = 1) where {T<:AbstractFloat}
    endo   = m.endogenous_states
    pseudo = m.pseudo_observables

    _n_states = n_states_augmented(m)
    _n_pseudo = n_pseudo_observables(m)

    # Initialize pseudo ZZ and DD matrices
    ZZ_pseudo = zeros(_n_pseudo, _n_states)
    DD_pseudo = zeros(_n_pseudo)

    ##########################################################
    ## PSEUDO-OBSERVABLE EQUATIONS
    ##########################################################

    ## Output
    ZZ_pseudo[pseudo[:y_t],endo[:y_t]] = 1.

    ## π_t
    ZZ_pseudo[pseudo[:π_t],endo[:π_t]] = 1.
    DD_pseudo[pseudo[:π_t]]            = 100*(m[:π_star]-1);

    ## z_t
    ZZ_pseudo[pseudo[:z_t], endo[:z_t]] = 1.

    ## Nominal FFR
    ZZ_pseudo[pseudo[:NominalFFR], endo[:R_t]] = 1.
    DD_pseudo[pseudo[:NominalFFR]] = m[:π_star] + m[:rA] + 4.0*m[:γ_Q]

    ## Real FFR
    ZZ_pseudo[pseudo[:RealFFR], endo[:R_t]] = 1.
    ZZ_pseudo[pseudo[:RealFFR], endo[:π_t]] = -4.
    DD_pseudo[pseudo[:RealFFR]] = m[:π_star] + m[:rA] + 4.0*m[:γ_Q] - 4.0*(100. * (m[:π_star] - 1.))

    return PseudoMeasurement(ZZ_pseudo, DD_pseudo)
end
