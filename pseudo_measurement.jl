 """
```
pseudo_measurement{T<:AbstractFloat}(m::Model1010{T})
```

Assign pseudo-measurement equation (a linear combination of states):

```
X_t = ZZ_pseudo*S_t + DD_pseudo
```
"""
function pseudo_measurement{T<:AbstractFloat}(m::Model1010{T})

    endo      = m.endogenous_states
    endo_addl = m.endogenous_states_augmented
    pseudo_names = [:y_t, :y_f_t, :NaturalRate, :π_t, :OutputGap, :ExAnteRealRate, :LongRunInflation,
                    :MarginalCost, :Wages, :FlexibleWages, :Hours, :FlexibleHours, :z_t,
                    :Expected10YearRateGap, :NominalFFR, :Expected10YearNominalRate,
                    :Expected10YearNominalNaturalRate,
                    :ExpectedNominalNaturalRate, :NominalRateGap, :LaborProductivityGrowth,
                    :Expected5YearRealRate, :Expected5YearRealNaturalRate, :Expected5YearRateGap,
                    :Expected5YearNominalRate, :Expected5YearNominalNaturalRate, :Expected10YearRealRate,
                    :Expected10YearRealNaturalRate]

    # Map pseudoobservables to indices
    pseudo_inds = Dict{Symbol,Int}()
    for (i,k) in enumerate(pseudo_names)
        pseudo_inds[k] = i
    end

    # Create PseudoObservable objects
    pseudo = Dict{Symbol,PseudoObservable}()
    for k in pseudo_names
        pseudo[k] = PseudoObservable(k)
    end

    # Compute TTT^10, used for Expected10YearRateGap, Expected10YearRate, and Expected10YearNaturalRate
    TTT, _, _ = solve(m)
    TTT10 = (1/40)*((UniformScaling(1.) - TTT)\(UniformScaling(1.) - TTT^40))
    TTT5 = (1/20)*((UniformScaling(1.) - TTT)\(UniformScaling(1.) - TTT^20))


    # Initialize pseudo ZZ and DD matrices
    ZZ_pseudo = zeros(length(pseudo), n_states_augmented(m))
    DD_pseudo = zeros(length(pseudo))

    ##########################################################
    ## PSEUDO-OBSERVABLE EQUATIONS
    ##########################################################

    ## Output
    ZZ_pseudo[pseudo_inds[:y_t],endo[:y_t]] = 1.
    pseudo[:y_t].name = "Output Growth"
    pseudo[:y_t].longname = "Output Growth Per Capita"

    ## Flexible Output
    ZZ_pseudo[pseudo_inds[:y_f_t],endo[:y_f_t]] = 1.
    pseudo[:y_f_t].name = "Flexible Output Growth"
    pseudo[:y_f_t].longname = "Output that would prevail in a flexible-price economy."

    ## Natural Rate
    ZZ_pseudo[pseudo_inds[:NaturalRate],endo[:r_f_t]] = 1.
    DD_pseudo[pseudo_inds[:NaturalRate]]              = 100.*(m[:rstar]-1.)

    pseudo[:NaturalRate].name = "Real Natural Rate"
    pseudo[:NaturalRate].longname = "The real interest rate that would prevail in a flexible-price economy."
    pseudo[:NaturalRate].rev_transform = quartertoannual

    ## π_t
    ZZ_pseudo[pseudo_inds[:π_t],endo[:π_t]] = 1.
    DD_pseudo[pseudo_inds[:π_t]]            = 100*(m[:π_star]-1);

    pseudo[:π_t].name = "Inflation"
    pseudo[:π_t].longname = "Inflation"
    pseudo[:π_t].rev_transform = quartertoannual

    ## Output Gap
    ZZ_pseudo[pseudo_inds[:OutputGap],endo[:y_t]] = 1;
    ZZ_pseudo[pseudo_inds[:OutputGap],endo[:y_f_t]] = -1;

    pseudo[:OutputGap].name = "Output Gap"
    pseudo[:OutputGap].longname = "Output Gap"

    ## Ex Ante Real Rate
    ZZ_pseudo[pseudo_inds[:ExAnteRealRate],endo[:R_t]]  = 1;
    ZZ_pseudo[pseudo_inds[:ExAnteRealRate],endo[:Eπ_t]] = -1;
    DD_pseudo[pseudo_inds[:ExAnteRealRate]]             = m[:Rstarn] - 100*(m[:π_star]-1);

    pseudo[:ExAnteRealRate].name = "Ex Ante Real Rate"
    pseudo[:ExAnteRealRate].longname = "Ex Ante Real Rate"
    pseudo[:ExAnteRealRate].rev_transform = quartertoannual

    ## Long Run Inflation
    ZZ_pseudo[pseudo_inds[:LongRunInflation],endo[:π_star_t]] = 1.
    DD_pseudo[pseudo_inds[:LongRunInflation]]                 = 100. *(m[:π_star]-1.)

    pseudo[:LongRunInflation].name = "Long Run Inflation"
    pseudo[:LongRunInflation].longname = "Long Run Inflation"
    pseudo[:LongRunInflation].rev_transform = quartertoannual

    ## Marginal Cost
    ZZ_pseudo[pseudo_inds[:MarginalCost],endo[:mc_t]] = 1.

    pseudo[:MarginalCost].name = "Marginal Cost"
    pseudo[:MarginalCost].longname = "Marginal Cost"

    ## Wages
    ZZ_pseudo[pseudo_inds[:Wages],endo[:w_t]] = 1.

    pseudo[:Wages].name = "Wages"
    pseudo[:Wages].longname = "Wages"

    ## Flexible Wages
    ZZ_pseudo[pseudo_inds[:FlexibleWages],endo[:w_f_t]] = 1.

    pseudo[:FlexibleWages].name = "Flexible Wages"
    pseudo[:FlexibleWages].longname = "Wages that would prevail in a flexible-wage economy"

    ## Hours
    ZZ_pseudo[pseudo_inds[:Hours],endo[:L_t]] = 1.

    pseudo[:Hours].name = "Hours"
    pseudo[:Hours].longname = "Hours"

    ## Flexible Hours
    ZZ_pseudo[pseudo_inds[:FlexibleHours],endo[:L_f_t]] = 1.
    pseudo[:FlexibleHours].name     = "Flexible Hours"
    pseudo[:FlexibleHours].longname = "Flexible Hours"

    ## z_t
    ZZ_pseudo[pseudo_inds[:z_t], endo[:z_t]] = 1.
    pseudo[:z_t].name     = "z_t"
    pseudo[:z_t].longname = "z_t"

    ## Expected 10-Year Rate Gap
    ZZ_pseudo[pseudo_inds[:Expected10YearRateGap], :] = TTT10[endo[:R_t], :] - TTT10[endo[:r_f_t], :] - TTT10[endo[:Eπ_t], :]
    pseudo[:Expected10YearRateGap].name     = "Expected 10-Year Rate Gap"
    pseudo[:Expected10YearRateGap].longname = "Expected 10-Year Rate Gap"
    pseudo[:Expected10YearRateGap].rev_transform = quartertoannual

    ## Nominal FFR
    ZZ_pseudo[pseudo_inds[:NominalFFR], endo[:R_t]] = 1.
    DD_pseudo[pseudo_inds[:NominalFFR]] = m[:Rstarn]
    pseudo[:NominalFFR].name     = "Nominal FFR"
    pseudo[:NominalFFR].longname = "Nominal FFR at an annual rate"
    pseudo[:NominalFFR].rev_transform = quartertoannual

    ## Expected 10-Year Nominal Interest Rate
    ZZ_pseudo[pseudo_inds[:Expected10YearNominalRate], :] = TTT10[endo[:R_t], :]
    DD_pseudo[pseudo_inds[:Expected10YearNominalRate]]    = m[:Rstarn]
    pseudo[:Expected10YearNominalRate].name     = "Expected 10-Year Nominal Rate"
    pseudo[:Expected10YearNominalRate].longname = "Expected 10-Year Nominal Interest Rate"
    pseudo[:Expected10YearNominalRate].rev_transform = quartertoannual

    ## Expected 10-Year Nominal Natural Rate
    ZZ_pseudo[pseudo_inds[:Expected10YearNominalNaturalRate], :] = TTT10[endo[:r_f_t], :] + TTT10[endo[:Eπ_t], :]
    DD_pseudo[pseudo_inds[:Expected10YearNominalNaturalRate]]    = m[:Rstarn]
    pseudo[:Expected10YearNominalNaturalRate].name     = "Expected 10-Year Nominal Natural Rate"
    pseudo[:Expected10YearNominalNaturalRate].longname = "Expected 10-Year Nominal Natural Rate of Interest"
    pseudo[:Expected10YearNominalNaturalRate].rev_transform = quartertoannual

    ## Expected Nominal Natural Rate
    ZZ_pseudo[pseudo_inds[:ExpectedNominalNaturalRate], endo[:r_f_t]] = 1.
    ZZ_pseudo[pseudo_inds[:ExpectedNominalNaturalRate], endo[:Eπ_t]]  = 1.
    DD_pseudo[pseudo_inds[:ExpectedNominalNaturalRate]]               = m[:Rstarn]
    pseudo[:ExpectedNominalNaturalRate].name     = "Expected Nominal Natural Rate"
    pseudo[:ExpectedNominalNaturalRate].longname = "Natural Rate + Expected Inflation"
    pseudo[:ExpectedNominalNaturalRate].rev_transform = quartertoannual

    ## Nominal Rate Gap
    ZZ_pseudo[pseudo_inds[:NominalRateGap], endo[:R_t]]   = 1.
    ZZ_pseudo[pseudo_inds[:NominalRateGap], endo[:r_f_t]] = -1.
    ZZ_pseudo[pseudo_inds[:NominalRateGap], endo[:Eπ_t]]  = -1.
    pseudo[:NominalRateGap].name     = "Nominal Rate Gap"
    pseudo[:NominalRateGap].longname = "Nominal FFR - Nominal Natural Rate"
    pseudo[:NominalRateGap].rev_transform = quartertoannual

    ## Labor Productivity Growth
    ZZ_pseudo[pseudo_inds[:LaborProductivityGrowth], endo[:y_t]]           = 1.
    ZZ_pseudo[pseudo_inds[:LaborProductivityGrowth], endo_addl[:y_t1]]     = -1.
    ZZ_pseudo[pseudo_inds[:LaborProductivityGrowth], endo[:z_t]]           = 1.
    ZZ_pseudo[pseudo_inds[:LaborProductivityGrowth], endo_addl[:e_gdp_t]]  = 1.
    ZZ_pseudo[pseudo_inds[:LaborProductivityGrowth], endo_addl[:e_gdp_t1]] = -m[:me_level]
    ZZ_pseudo[pseudo_inds[:LaborProductivityGrowth], endo[:L_t]]           = -1
    ZZ_pseudo[pseudo_inds[:LaborProductivityGrowth], endo_addl[:L_t1]]     = 1.
    DD_pseudo[pseudo_inds[:LaborProductivityGrowth]]                       = 100*(exp(m[:z_star]) - 1)
    pseudo[:LaborProductivityGrowth].name     = "Labor Productivity Growth"
    pseudo[:LaborProductivityGrowth].longname = "Labor Productivity Growth Rate"
    pseudo[:LaborProductivityGrowth].rev_transform = quartertoannual

    ## Expected 5-Year Real Interest Rate
    ZZ_pseudo[pseudo_inds[:Expected5YearRealRate], :] = TTT5[endo[:R_t], :] - TTT5[endo[:Eπ_t], :]
    DD_pseudo[pseudo_inds[:Expected5YearRealRate]]    = m[:Rstarn] - 100*(m[:π_star]-1);

    pseudo[:Expected5YearRealRate].name     = "Expected 5-Year Real Rate"
    pseudo[:Expected5YearRealRate].longname = "Expected 5-Year Real Interest Rate"
    pseudo[:Expected5YearRealRate].rev_transform = quartertoannual

    ## Expected 5-Year Real Natural Rate
    ZZ_pseudo[pseudo_inds[:Expected5YearRealNaturalRate], :] = TTT5[endo[:r_f_t], :]
    DD_pseudo[pseudo_inds[:Expected5YearRealNaturalRate]]    = m[:Rstarn] - 100*(m[:π_star]-1);
    pseudo[:Expected5YearRealNaturalRate].name     = "Expected 5-Year Real Natural Rate"
    pseudo[:Expected5YearRealNaturalRate].longname = "Expected 5-Year Real Natural Rate of Interest"
    pseudo[:Expected5YearRealNaturalRate].rev_transform = quartertoannual

    ## Expected 5-Year Rate Gap
    ZZ_pseudo[pseudo_inds[:Expected5YearRateGap], :] = TTT5[endo[:R_t], :] -
        TTT5[endo[:r_f_t], :] - TTT5[endo[:Eπ_t], :]
    pseudo[:Expected5YearRateGap].name     = "Expected 5-Year Rate Gap"
    pseudo[:Expected5YearRateGap].longname = "Expected 5-Year Rate Gap"
    pseudo[:Expected5YearRateGap].rev_transform = quartertoannual

    ## Expected 5-Year Nominal Interest Rate
    ZZ_pseudo[pseudo_inds[:Expected5YearNominalRate], :] = TTT5[endo[:R_t], :]
    DD_pseudo[pseudo_inds[:Expected5YearNominalRate]]    = m[:Rstarn]
    pseudo[:Expected5YearNominalRate].name     = "Expected 5-Year Rate"
    pseudo[:Expected5YearNominalRate].longname = "Expected 5-Year Interest Rate"
    pseudo[:Expected5YearNominalRate].rev_transform = quartertoannual

    ## Expected 5-Year Nominal Natural Rate
    ZZ_pseudo[pseudo_inds[:Expected5YearNominalNaturalRate], :] = TTT5[endo[:r_f_t], :] + TTT5[endo[:Eπ_t], :]
    DD_pseudo[pseudo_inds[:Expected5YearNominalNaturalRate]]    = m[:Rstarn]
    pseudo[:Expected5YearNominalNaturalRate].name     = "Expected 5-Year Natural Rate"
    pseudo[:Expected5YearNominalNaturalRate].longname = "Expected 5-Year Natural Rate of Interest"
    pseudo[:Expected5YearNominalNaturalRate].rev_transform = quartertoannual

    ## Expected 10-Year Real Interest Rate
    ZZ_pseudo[pseudo_inds[:Expected10YearRealRate], :] = TTT10[endo[:R_t], :] - TTT10[endo[:Eπ_t], :]
    DD_pseudo[pseudo_inds[:Expected10YearRealRate]]    = m[:Rstarn] - 100*(m[:π_star]-1);

    pseudo[:Expected10YearRealRate].name     = "Expected 10-Year Real Rate"
    pseudo[:Expected10YearRealRate].longname = "Expected 10-Year Real Interest Rate"
    pseudo[:Expected10YearRealRate].rev_transform = quartertoannual

    ## Expected 10-Year Real Natural Rate
    ZZ_pseudo[pseudo_inds[:Expected10YearRealNaturalRate], :] = TTT10[endo[:r_f_t], :]
    DD_pseudo[pseudo_inds[:Expected10YearRealNaturalRate]]    = m[:Rstarn] - 100*(m[:π_star]-1);
    pseudo[:Expected10YearRealNaturalRate].name     = "Expected 10-Year Real Natural Rate"
    pseudo[:Expected10YearRealNaturalRate].longname = "Expected 10-Year Real Natural Rate of Interest"
    pseudo[:Expected10YearRealNaturalRate].rev_transform = quartertoannual

    # Collect indices and transforms
    return pseudo, PseudoObservableMapping(pseudo_inds, ZZ_pseudo, DD_pseudo)
end