"""
```
pseudo_measurement{T<:AbstractFloat}(m::Model990{T})
```

Assign pseudo-measurement equation (a linear combination of states):

```
X_t = ZZ_pseudo*S_t + DD_pseudo
```
"""
function pseudo_measurement{T<:AbstractFloat}(m::Model990{T}; apply_altpolicy::Bool = false)

    endo = m.endogenous_states
    pseudo_names = [:y_t, :y_f_t, :NaturalRate, :π_t, :OutputGap, :ExAnteRealRate, :LongRunInflation,
                    :MarginalCost, :Wages, :FlexibleWages, :Hours, :FlexibleHours, :z_t, :Z_t, :ztil_t,
                    :NominalFFR, :RealFFR, :ExpectedNominalNaturalRate]

    # Map pseudoobservables to indices
    pseudo_inds = OrderedDict{Symbol,Int}()
    for (i,k) in enumerate(pseudo_names)
        pseudo_inds[k] = i
    end

    # Create PseudoObservable objects
    pseudo = OrderedDict{Symbol,PseudoObservable}()
    for k in pseudo_names
        pseudo[k] = PseudoObservable(k)
    end

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
    pseudo[:y_f_t].longname = "Output that would obtain in a flexible-price economy."

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
    pseudo[:FlexibleWages].longname = "Wages that would obtain in a flexible-wage economy"

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

    ## Z_t
    ZZ_pseudo[pseudo_inds[:Z_t], endo[:z_t]] = 1.
    DD_pseudo[pseudo_inds[:Z_t]]             = 100. * exp(m[:z_star]-1)
    pseudo[:Z_t].name     = "Z_t"
    pseudo[:Z_t].longname = "Z_t"

    ## ztil_t
    ZZ_pseudo[pseudo_inds[:ztil_t], endo[:ztil_t]] = 1.
    pseudo[:z_t].name     = "ztil_t"
    pseudo[:z_t].longname = "ztil_t"

    ## Nominal FFR
    ZZ_pseudo[pseudo_inds[:NominalFFR], endo[:R_t]] = 1.
    DD_pseudo[pseudo_inds[:NominalFFR]] = m[:Rstarn]
    pseudo[:NominalFFR].name     = "Nominal FFR"
    pseudo[:NominalFFR].longname = "Nominal FFR at an annual rate"
    pseudo[:NominalFFR].rev_transform = quartertoannual

    ## Real FFR
    ZZ_pseudo[pseudo_inds[:RealFFR], endo[:R_t]] = 1.
    ZZ_pseudo[pseudo_inds[:RealFFR], endo[:π_t]] = -1.
    DD_pseudo[pseudo_inds[:RealFFR]] =  m[:Rstarn] - (100. * (m[:π_star] - 1.))
    pseudo[:RealFFR].name     = "Real FFR"
    pseudo[:RealFFR].longname = "Real FFR at an annual rate"
    pseudo[:RealFFR].rev_transform = quartertoannual

    ## Expected nominal natural rate
    ZZ_pseudo[pseudo_inds[:ExpectedNominalNaturalRate], endo[:r_f_t]] = 1.
    ZZ_pseudo[pseudo_inds[:ExpectedNominalNaturalRate], endo[:Eπ_t]]  = 1.
    DD_pseudo[pseudo_inds[:ExpectedNominalNaturalRate]]               =  m[:Rstarn]
    pseudo[:ExpectedNominalNaturalRate].name  = "Expected Nominal Natural Rate"
    pseudo[:ExpectedNominalNaturalRate].longname = "Natural Rate + Expected Inflation"
    pseudo[:ExpectedNominalNaturalRate].rev_transform = quartertoannual

    # Collect indices and transforms
    return pseudo, PseudoObservableMapping(pseudo_inds, ZZ_pseudo, DD_pseudo)
end

