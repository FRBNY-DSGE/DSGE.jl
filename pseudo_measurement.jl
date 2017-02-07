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
                    :ExpectedAvg10YearRateGap, :NominalFFR, :ExpectedAvg10YearNominalRate,
                    :ExpectedAvg10YearNominalNaturalRate,
                    :ExpectedAvgNominalNaturalRate, :NominalRateGap, :LaborProductivityGrowth,
                    :ExpectedAvg5YearRealRate,     :ExpectedAvg5YearRealNaturalRate, :ExpectedAvg5YearRateGap,
                    :ExpectedAvg5YearNominalRate,  :ExpectedAvg5YearNominalNaturalRate,
                    :ExpectedAvg10YearRealRate,    :ExpectedAvg10YearRealNaturalRate,
                    :ExpectedAvg20YearNominalRate, :ExpectedAvg20YearNominalNaturalRate,
                    :ExpectedAvg20YearRealRate,    :ExpectedAvg20YearRealNaturalRate,
                    :ExpectedAvg20YearRateGap,
                    # :Forward5YearRealNaturalRate,
                    # :Forward10YearRealNaturalRate,
                    :Forward20YearRealNaturalRate,
                    :Forward25YearRealNaturalRate,
                    :Forward30YearRealNaturalRate]

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

    # Compute TTT^10, used for ExpectedAvg10YearRateGap, ExpectedAvg10YearRate, and ExpectedAvg10YearNaturalRate
    TTT, _, _ = solve(m)
    TTT20 = (1/80)*((UniformScaling(1.) - TTT)\(UniformScaling(1.) - TTT^80))
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

    ## Expected Average 10-Year Rate Gap
    ZZ_pseudo[pseudo_inds[:ExpectedAvg10YearRateGap], :] = TTT10[endo[:R_t], :] - TTT10[endo[:r_f_t], :] - TTT10[endo[:Eπ_t], :]
    pseudo[:ExpectedAvg10YearRateGap].name     = "Expected Average 10-Year Rate Gap"
    pseudo[:ExpectedAvg10YearRateGap].longname = "Expected Average 10-Year Rate Gap"
    pseudo[:ExpectedAvg10YearRateGap].rev_transform = quartertoannual

    ## Nominal FFR
    ZZ_pseudo[pseudo_inds[:NominalFFR], endo[:R_t]] = 1.
    DD_pseudo[pseudo_inds[:NominalFFR]] = m[:Rstarn]
    pseudo[:NominalFFR].name     = "Nominal FFR"
    pseudo[:NominalFFR].longname = "Nominal FFR at an annual rate"
    pseudo[:NominalFFR].rev_transform = quartertoannual

    ## Expected Average 10-Year Nominal Interest Rate
    ZZ_pseudo[pseudo_inds[:ExpectedAvg10YearNominalRate], :] = TTT10[endo[:R_t], :]
    DD_pseudo[pseudo_inds[:ExpectedAvg10YearNominalRate]]    = m[:Rstarn]
    pseudo[:ExpectedAvg10YearNominalRate].name     = "Expected Average 10-Year Nominal Rate"
    pseudo[:ExpectedAvg10YearNominalRate].longname = "Expected Average 10-Year Nominal Interest Rate"
    pseudo[:ExpectedAvg10YearNominalRate].rev_transform = quartertoannual

    ## Expected Average 10-Year Nominal Natural Rate
    ZZ_pseudo[pseudo_inds[:ExpectedAvg10YearNominalNaturalRate], :] = TTT10[endo[:r_f_t], :] + TTT10[endo[:Eπ_t], :]
    DD_pseudo[pseudo_inds[:ExpectedAvg10YearNominalNaturalRate]]    = m[:Rstarn]
    pseudo[:ExpectedAvg10YearNominalNaturalRate].name     = "Expected Average 10-Year Nominal Natural Rate"
    pseudo[:ExpectedAvg10YearNominalNaturalRate].longname = "Expected Average 10-Year Nominal Natural Rate of Interest"
    pseudo[:ExpectedAvg10YearNominalNaturalRate].rev_transform = quartertoannual

    ## Expected AverageNominal Natural Rate
    ZZ_pseudo[pseudo_inds[:ExpectedAvgNominalNaturalRate], endo[:r_f_t]] = 1.
    ZZ_pseudo[pseudo_inds[:ExpectedAvgNominalNaturalRate], endo[:Eπ_t]]  = 1.
    DD_pseudo[pseudo_inds[:ExpectedAvgNominalNaturalRate]]               = m[:Rstarn]
    pseudo[:ExpectedAvgNominalNaturalRate].name     = "Expected Average Nominal Natural Rate"
    pseudo[:ExpectedAvgNominalNaturalRate].longname = "Natural Rate + Expected Inflation"
    pseudo[:ExpectedAvgNominalNaturalRate].rev_transform = quartertoannual

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

    ## Expected Average 5-Year Real Interest Rate
    ZZ_pseudo[pseudo_inds[:ExpectedAvg5YearRealRate], :] = TTT5[endo[:R_t], :] - TTT5[endo[:Eπ_t], :]
    DD_pseudo[pseudo_inds[:ExpectedAvg5YearRealRate]]    = m[:Rstarn] - 100*(m[:π_star]-1);

    pseudo[:ExpectedAvg5YearRealRate].name     = "Expected Average 5-Year Real Rate"
    pseudo[:ExpectedAvg5YearRealRate].longname = "Expected Average 5-Year Real Interest Rate"
    pseudo[:ExpectedAvg5YearRealRate].rev_transform = quartertoannual

    ## Expected Average 5-Year Real Natural Rate
    ZZ_pseudo[pseudo_inds[:ExpectedAvg5YearRealNaturalRate], :] = TTT5[endo[:r_f_t], :]
    DD_pseudo[pseudo_inds[:ExpectedAvg5YearRealNaturalRate]]    = m[:Rstarn] - 100*(m[:π_star]-1);
    pseudo[:ExpectedAvg5YearRealNaturalRate].name     = "Expected Average 5-Year Real Natural Rate"
    pseudo[:ExpectedAvg5YearRealNaturalRate].longname = "Expected Average 5-Year Real Natural Rate of Interest"
    pseudo[:ExpectedAvg5YearRealNaturalRate].rev_transform = quartertoannual

    ## Expected Average 5-Year Rate Gap
    ZZ_pseudo[pseudo_inds[:ExpectedAvg5YearRateGap], :] = TTT5[endo[:R_t], :] -
        TTT5[endo[:r_f_t], :] - TTT5[endo[:Eπ_t], :]
    pseudo[:ExpectedAvg5YearRateGap].name     = "Expected Average 5-Year Rate Gap"
    pseudo[:ExpectedAvg5YearRateGap].longname = "Expected Average 5-Year Rate Gap"
    pseudo[:ExpectedAvg5YearRateGap].rev_transform = quartertoannual

    ## Expected Average 5-Year Nominal Interest Rate
    ZZ_pseudo[pseudo_inds[:ExpectedAvg5YearNominalRate], :] = TTT5[endo[:R_t], :]
    DD_pseudo[pseudo_inds[:ExpectedAvg5YearNominalRate]]    = m[:Rstarn]
    pseudo[:ExpectedAvg5YearNominalRate].name     = "Expected Average 5-Year Rate"
    pseudo[:ExpectedAvg5YearNominalRate].longname = "Expected Average 5-Year Interest Rate"
    pseudo[:ExpectedAvg5YearNominalRate].rev_transform = quartertoannual

    ## Expected Average 5-Year Nominal Natural Rate
    ZZ_pseudo[pseudo_inds[:ExpectedAvg5YearNominalNaturalRate], :] = TTT5[endo[:r_f_t], :] + TTT5[endo[:Eπ_t], :]
    DD_pseudo[pseudo_inds[:ExpectedAvg5YearNominalNaturalRate]]    = m[:Rstarn]
    pseudo[:ExpectedAvg5YearNominalNaturalRate].name     = "Expected Average 5-Year Natural Rate"
    pseudo[:ExpectedAvg5YearNominalNaturalRate].longname = "Expected Average 5-Year Natural Rate of Interest"
    pseudo[:ExpectedAvg5YearNominalNaturalRate].rev_transform = quartertoannual

    ## Expected Average 10-Year Real Interest Rate
    ZZ_pseudo[pseudo_inds[:ExpectedAvg10YearRealRate], :] = TTT10[endo[:R_t], :] - TTT10[endo[:Eπ_t], :]
    DD_pseudo[pseudo_inds[:ExpectedAvg10YearRealRate]]    = m[:Rstarn] - 100*(m[:π_star]-1);

    pseudo[:ExpectedAvg10YearRealRate].name     = "Expected Average 10-Year Real Rate"
    pseudo[:ExpectedAvg10YearRealRate].longname = "Expected Average 10-Year Real Interest Rate"
    pseudo[:ExpectedAvg10YearRealRate].rev_transform = quartertoannual

    ## Expected Average 10-Year Real Natural Rate
    ZZ_pseudo[pseudo_inds[:ExpectedAvg10YearRealNaturalRate], :] = TTT10[endo[:r_f_t], :]
    DD_pseudo[pseudo_inds[:ExpectedAvg10YearRealNaturalRate]]    = m[:Rstarn] - 100*(m[:π_star]-1);
    pseudo[:ExpectedAvg10YearRealNaturalRate].name     = "Expected Average 10-Year Real Natural Rate"
    pseudo[:ExpectedAvg10YearRealNaturalRate].longname = "Expected Average 10-Year Real Natural Rate of Interest"
    pseudo[:ExpectedAvg10YearRealNaturalRate].rev_transform = quartertoannual

    ## Expected Average 20-Year Nominal Interest Rate
    ZZ_pseudo[pseudo_inds[:ExpectedAvg20YearNominalRate], :] = TTT20[endo[:R_t], :]
    DD_pseudo[pseudo_inds[:ExpectedAvg20YearNominalRate]]    = m[:Rstarn]
    pseudo[:ExpectedAvg20YearNominalRate].name     = "Expected Average 20-Year Nominal Rate"
    pseudo[:ExpectedAvg20YearNominalRate].longname = "Expected Average 20-Year Nominal Interest Rate"
    pseudo[:ExpectedAvg20YearNominalRate].rev_transform = quartertoannual

    ## Expected Average 20-Year Nominal Natural Rate
    ZZ_pseudo[pseudo_inds[:ExpectedAvg20YearNominalNaturalRate], :] = TTT20[endo[:r_f_t], :] + TTT20[endo[:Eπ_t], :]
    DD_pseudo[pseudo_inds[:ExpectedAvg20YearNominalNaturalRate]]    = m[:Rstarn]
    pseudo[:ExpectedAvg20YearNominalNaturalRate].name     = "Expected Average 20-Year Nominal Natural Rate"
    pseudo[:ExpectedAvg20YearNominalNaturalRate].longname = "Expected Average 20-Year Nominal Natural Rate of Interest"
    pseudo[:ExpectedAvg20YearNominalNaturalRate].rev_transform = quartertoannual

    ## Expected Average 20-Year Real Interest Rate
    ZZ_pseudo[pseudo_inds[:ExpectedAvg20YearRealRate], :] = TTT20[endo[:R_t], :] - TTT20[endo[:Eπ_t], :]
    DD_pseudo[pseudo_inds[:ExpectedAvg20YearRealRate]]    = m[:Rstarn] - 100*(m[:π_star]-1);

    pseudo[:ExpectedAvg20YearRealRate].name     = "Expected Average 20-Year Real Rate"
    pseudo[:ExpectedAvg20YearRealRate].longname = "Expected Average 20-Year Real Interest Rate"
    pseudo[:ExpectedAvg20YearRealRate].rev_transform = quartertoannual

    ## Expected Average 20-Year Real Natural Rate
    ZZ_pseudo[pseudo_inds[:ExpectedAvg20YearRealNaturalRate], :] = TTT20[endo[:r_f_t], :]
    DD_pseudo[pseudo_inds[:ExpectedAvg20YearRealNaturalRate]]    = m[:Rstarn] - 100*(m[:π_star]-1);
    pseudo[:ExpectedAvg20YearRealNaturalRate].name     = "Expected Average 20-Year Real Natural Rate"
    pseudo[:ExpectedAvg20YearRealNaturalRate].longname = "Expected Average 20-Year Real Natural Rate of Interest"
    pseudo[:ExpectedAvg20YearRealNaturalRate].rev_transform = quartertoannual

    ## Expected Average 20-Year Rate Gap
    ZZ_pseudo[pseudo_inds[:ExpectedAvg20YearRateGap], :] = TTT20[endo[:R_t], :] -
        TTT20[endo[:r_f_t], :] - TTT20[endo[:Eπ_t], :]
    pseudo[:ExpectedAvg20YearRateGap].name     = "Expected Average 20-Year Rate Gap"
    pseudo[:ExpectedAvg20YearRateGap].longname = "Expected Average 20-Year Rate Gap"
    pseudo[:ExpectedAvg20YearRateGap].rev_transform = quartertoannual

    # ## 5-Year Real Natural Forward Rate
    # TTT5_sim = TTT^20
    # ZZ_pseudo[pseudo_inds[:Forward5YearRealNaturalRate], :] = ZZ_pseudo[pseudo_inds[:NaturalRate], :] * TTT5_sim
    # DD_pseudo[pseudo_inds[:Forward5YearRealNaturalRate]]    = m[:Rstarn] - 100*(m[:π_star]-1);
    # pseudo[:Forward5YearRealNaturalRate].name     = "Forward 5-Year Real Natural Rate"
    # pseudo[:Forward5YearRealNaturalRate].longname = "Forward 5-Year Real Natural Rate of Interest (not the average, computed by projecting TTT foreward 80 periods)"
    # pseudo[:Forward5YearRealNaturalRate].rev_transform = quartertoannual

    # ## 10-Year Real Natural Forward Rate
    # TTT10_sim = TTT^40
    # ZZ_pseudo[pseudo_inds[:Forward10YearRealNaturalRate], :] = ZZ_pseudo[pseudo_inds[:NaturalRate], :] * TTT10_sim
    # DD_pseudo[pseudo_inds[:Forward10YearRealNaturalRate]]    = m[:Rstarn] - 100*(m[:π_star]-1);
    # pseudo[:Forward10YearRealNaturalRate].name     = "Forward 10-Year Real Natural Rate"
    # pseudo[:Forward10YearRealNaturalRate].longname = "Forward 10-Year Real Natural Rate of Interest (not the average, computed by projecting TTT foreward 80 periods)"
    # pseudo[:Forward10YearRealNaturalRate].rev_transform = quartertoannual

    ## 20-Year Real Natural Forward Rate
    TTT20_sim = TTT^80
    ZZ_pseudo[pseudo_inds[:Forward20YearRealNaturalRate], :] = ZZ_pseudo[pseudo_inds[:NaturalRate], :] * TTT20_sim
    DD_pseudo[pseudo_inds[:Forward20YearRealNaturalRate]]    = m[:Rstarn] - 100*(m[:π_star]-1);
    pseudo[:Forward20YearRealNaturalRate].name     = "Forward 20-Year Real Natural Rate"
    pseudo[:Forward20YearRealNaturalRate].longname = "Forward 20-Year Real Natural Rate of Interest (not the average, computed by projecting TTT foreward 80 periods)"
    pseudo[:Forward20YearRealNaturalRate].rev_transform = quartertoannual

    ## 25-Year Real Natural Forward Rate
    TTT25_sim = TTT^100
    ZZ_pseudo[pseudo_inds[:Forward25YearRealNaturalRate], :] = ZZ_pseudo[pseudo_inds[:NaturalRate], :] * TTT25_sim
    DD_pseudo[pseudo_inds[:Forward25YearRealNaturalRate]]    = m[:Rstarn] - 100*(m[:π_star]-1);
    pseudo[:Forward25YearRealNaturalRate].name     = "Forward 25-Year Real Natural Rate"
    pseudo[:Forward25YearRealNaturalRate].longname = "Forward 25-Year Real Natural Rate of Interest (not the average, computed by projecting TTT foreward 80 periods)"
    pseudo[:Forward25YearRealNaturalRate].rev_transform = quartertoannual

    ## 30-Year Real Natural Forward Rate
    TTT30_sim = TTT^120
    ZZ_pseudo[pseudo_inds[:Forward30YearRealNaturalRate], :] = ZZ_pseudo[pseudo_inds[:NaturalRate], :] * TTT30_sim
    DD_pseudo[pseudo_inds[:Forward30YearRealNaturalRate]]    = m[:Rstarn] - 100*(m[:π_star]-1);
    pseudo[:Forward30YearRealNaturalRate].name     = "Forward 30-Year Real Natural Rate"
    pseudo[:Forward30YearRealNaturalRate].longname = "Forward 30-Year Real Natural Rate of Interest (not the average, computed by projecting TTT foreward 80 periods)"
    pseudo[:Forward30YearRealNaturalRate].rev_transform = quartertoannual


    # Collect indices and transforms
    return pseudo, PseudoObservableMapping(pseudo_inds, ZZ_pseudo, DD_pseudo)
end