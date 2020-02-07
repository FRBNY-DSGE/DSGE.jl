"""
```
pseudo_measurement{T<:AbstractFloat}(m::Model904{T},
    TTT::Matrix{T}, RRR::Matrix{T}, CCC::Vector{T})
```

Assign pseudo-measurement equation (a linear combination of states):

```
x_t = ZZ_pseudo*s_t + DD_pseudo
```
"""
function pseudo_measurement(m::Model904{T},
                                              TTT::Matrix{T},
                                              RRR::Matrix{T},
                                              CCC::Vector{T}) where {T<:AbstractFloat}
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

    ## Flexible Output
    ZZ_pseudo[pseudo[:y_f_t],endo[:y_f_t]] = 1.

    ## Natural Rate
    ZZ_pseudo[pseudo[:NaturalRate],endo[:r_f_t]] = 1.
    DD_pseudo[pseudo[:NaturalRate]]              = 100 .* (m[:rstar]-1.)

    ## π_t
    ZZ_pseudo[pseudo[:π_t],endo[:π_t]] = 1.
    DD_pseudo[pseudo[:π_t]]            = 100*(m[:π_star]-1);

    ## Output Gap
    ZZ_pseudo[pseudo[:OutputGap],endo[:y_t]] = 1;
    ZZ_pseudo[pseudo[:OutputGap],endo[:y_f_t]] = -1;

    ## Ex Ante Real Rate
    ZZ_pseudo[pseudo[:ExAnteRealRate],endo[:R_t]]  = 1;
    ZZ_pseudo[pseudo[:ExAnteRealRate],endo[:Eπ_t]] = -1;
    DD_pseudo[pseudo[:ExAnteRealRate]]             = m[:Rstarn] - 100*(m[:π_star]-1);

    ## Long Run Inflation
    ZZ_pseudo[pseudo[:LongRunInflation],endo[:π_star_t]] = 1.
    DD_pseudo[pseudo[:LongRunInflation]]                 = 100. *(m[:π_star]-1.)

    ## Marginal Cost
    ZZ_pseudo[pseudo[:MarginalCost],endo[:mc_t]] = 1.

    ## Wages
    ZZ_pseudo[pseudo[:Wages],endo[:w_t]] = 1.

    ## Flexible Wages
    ZZ_pseudo[pseudo[:FlexibleWages],endo[:w_f_t]] = 1.

    ## Hours
    ZZ_pseudo[pseudo[:Hours],endo[:L_t]] = 1.

    ## Flexible Hours
    ZZ_pseudo[pseudo[:FlexibleHours],endo[:L_f_t]] = 1.

    ## z_t
    ZZ_pseudo[pseudo[:z_t], endo[:z_t]] = 1.

    ## Z_t
    ZZ_pseudo[pseudo[:Z_t], endo[:z_t]] = 1.
    DD_pseudo[pseudo[:Z_t]]             = 100. * exp(m[:z_star]-1)

    ## ztil_t
    ZZ_pseudo[pseudo[:ztil_t], endo[:ztil_t]] = 1.

    ## Nominal FFR
    ZZ_pseudo[pseudo[:NominalFFR], endo[:R_t]] = 1.
    DD_pseudo[pseudo[:NominalFFR]] = m[:Rstarn]

    ## Real FFR
    ZZ_pseudo[pseudo[:RealFFR], endo[:R_t]] = 1.
    ZZ_pseudo[pseudo[:RealFFR], endo[:π_t]] = -1.
    DD_pseudo[pseudo[:RealFFR]] =  m[:Rstarn] - (100. * (m[:π_star] - 1.))

    ## Expected nominal natural rate
    ZZ_pseudo[pseudo[:ExpectedNominalNaturalRate], endo[:r_f_t]] = 1.
    ZZ_pseudo[pseudo[:ExpectedNominalNaturalRate], endo[:Eπ_t]]  = 1.
    DD_pseudo[pseudo[:ExpectedNominalNaturalRate]]               =  m[:Rstarn]

    ## S^Inf
    betabar     = (exp((1-m[:σ_c]) * m[:z_star]) * m[:β])
    S1 = inv(Matrix(1.0I, size(TTT, 1), size(TTT, 1)) - betabar*TTT)
    ZZ_pseudo[pseudo[:SInf], :] = S1[endo[:mc_t],:]
    ZZ_pseudo[pseudo[:SInf_t1], pseudo[:SInf]] = 1.0
    ZZ_pseudo[pseudo[:SInf_t2], pseudo[:SInf]] = 1.0
    ZZ_pseudo[pseudo[:SInf_t3], pseudo[:SInf]] = 1.0

    # π_til
    kappa =  ((1 - m[:ζ_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))*(1 - m[:ζ_p]))/(m[:ζ_p]*((m[:Φ]- 1)*m[:ϵ_p] + 1))/(1 + m[:ι_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
    exp( (m[:σ_c] -1) * m[:z_star]) / m[:β]
    kappa_coef = kappa * (1+m[:ι_p]*betabar)
    ZZ_pseudo[pseudo[:Pi_Til], pseudo[:SInf]] = kappa_coef
    ZZ_pseudo[pseudo[:Pi_Til], pseudo[:SInf_t1]] = kappa_coef*m[:ι_p]
    ZZ_pseudo[pseudo[:Pi_Til], pseudo[:SInf_t2]] = kappa_coef*m[:ι_p]^2
    ZZ_pseudo[pseudo[:Pi_Til], pseudo[:SInf_t3]] = kappa_coef*m[:ι_p]^3
    DD_pseudo[pseudo[:Pi_Til]] = 100*(m[:π_star] -1)

    # Collect indices and transforms
    return PseudoMeasurement(ZZ_pseudo, DD_pseudo)
end
