 """
```
pseudo_measurement{T<:AbstractFloat}(m::Model1002{T},
    TTT::Matrix{T}, RRR::Matrix{T}, CCC::Vector{T})
```

Assign pseudo-measurement equation (a linear combination of states):

```
x_t = ZZ_pseudo*s_t + DD_pseudo
```
"""
function pseudo_measurement{T<:AbstractFloat}(m::Model1002{T},
                                              TTT::Matrix{T},
                                              RRR::Matrix{T},
                                              CCC::Vector{T})

    endo      = m.endogenous_states
    endo_addl = m.endogenous_states_augmented
    pseudo    = m.pseudo_observables

    _n_states = n_states_augmented(m)
    _n_pseudo = n_pseudo_observables(m)

    # Compute TTT^10, used for Expected10YearRateGap, Expected10YearRate, and Expected10YearNaturalRate
    TTT10 = (1/40)*((UniformScaling(1.) - TTT)\(TTT - TTT^41))

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
    DD_pseudo[pseudo[:NaturalRate]]              = 100.*(m[:rstar]-1.)

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

    ## Expected 10-Year Rate Gap
    ZZ_pseudo[pseudo[:Expected10YearRateGap], :] = TTT10[endo[:R_t], :] - TTT10[endo[:r_f_t], :] - TTT10[endo[:Eπ_t], :]

    ## Nominal FFR
    ZZ_pseudo[pseudo[:NominalFFR], endo[:R_t]] = 1.
    DD_pseudo[pseudo[:NominalFFR]] = m[:Rstarn]

    ## Expected 10-Year Interest Rate
    ZZ_pseudo[pseudo[:Expected10YearRate], :] = TTT10[endo[:R_t], :]
    DD_pseudo[pseudo[:Expected10YearRate]]    = m[:Rstarn]

    ## Expected 10-Year Natural Rate
    ZZ_pseudo[pseudo[:Expected10YearNaturalRate], :] = TTT10[endo[:r_f_t], :] + TTT10[endo[:Eπ_t], :]
    DD_pseudo[pseudo[:Expected10YearNaturalRate]]    = m[:Rstarn]

    ## Expected Nominal Natural Rate
    ZZ_pseudo[pseudo[:ExpectedNominalNaturalRate], endo[:r_f_t]] = 1.
    ZZ_pseudo[pseudo[:ExpectedNominalNaturalRate], endo[:Eπ_t]]  = 1.
    DD_pseudo[pseudo[:ExpectedNominalNaturalRate]]               = m[:Rstarn]

    ## Nominal Rate Gap
    ZZ_pseudo[pseudo[:NominalRateGap], endo[:R_t]]   = 1.
    ZZ_pseudo[pseudo[:NominalRateGap], endo[:r_f_t]] = -1.
    ZZ_pseudo[pseudo[:NominalRateGap], endo[:Eπ_t]]  = -1.

    ## Labor Productivity Growth
    ZZ_pseudo[pseudo[:LaborProductivityGrowth], endo[:y_t]]           = 1.
    ZZ_pseudo[pseudo[:LaborProductivityGrowth], endo_addl[:y_t1]]     = -1.
    ZZ_pseudo[pseudo[:LaborProductivityGrowth], endo[:z_t]]           = 1.
    ZZ_pseudo[pseudo[:LaborProductivityGrowth], endo_addl[:e_gdp_t]]  = 1.
    ZZ_pseudo[pseudo[:LaborProductivityGrowth], endo_addl[:e_gdp_t1]] = -m[:me_level]
    ZZ_pseudo[pseudo[:LaborProductivityGrowth], endo[:L_t]]           = -1
    ZZ_pseudo[pseudo[:LaborProductivityGrowth], endo_addl[:L_t1]]     = 1.
    DD_pseudo[pseudo[:LaborProductivityGrowth]]                       = 100*(exp(m[:z_star]) - 1)

    return PseudoMeasurement(ZZ_pseudo, DD_pseudo)
end