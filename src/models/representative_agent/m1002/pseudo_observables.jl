function init_pseudo_observable_mappings!(m::Model1002)

    pseudo_names = [:y_t, :y_f_t, :NaturalRate, :π_t, :OutputGap, :ExAnteRealRate, :LongRunInflation,
                    :MarginalCost, :Wages, :FlexibleWages, :Hours, :FlexibleHours, :z_t,
                    :Expected10YearRateGap, :NominalFFR, :Expected10YearRate, :Expected10YearNaturalRate,
                    :ExpectedNominalNaturalRate, :NominalRateGap, :LaborProductivityGrowth]

    # Create PseudoObservable objects
    pseudo = OrderedDict{Symbol,PseudoObservable}()
    for k in pseudo_names
        pseudo[k] = PseudoObservable(k)
    end

    # Fill in names and reverse transforms
    pseudo[:y_t].name = "Output Growth"
    pseudo[:y_t].longname = "Output Growth Per Capita"

    pseudo[:y_f_t].name = "Flexible Output Growth"
    pseudo[:y_f_t].longname = "Output that would prevail in a flexible-price economy."

    pseudo[:NaturalRate].name = "Real Natural Rate"
    pseudo[:NaturalRate].longname = "The real interest rate that would prevail in a flexible-price economy."
    pseudo[:NaturalRate].rev_transform = quartertoannual

    pseudo[:π_t].name = "Inflation"
    pseudo[:π_t].longname = "Inflation"
    pseudo[:π_t].rev_transform = quartertoannual

    pseudo[:OutputGap].name = "Output Gap"
    pseudo[:OutputGap].longname = "Output Gap"

    pseudo[:ExAnteRealRate].name = "Ex Ante Real Rate"
    pseudo[:ExAnteRealRate].longname = "Ex Ante Real Rate"
    pseudo[:ExAnteRealRate].rev_transform = quartertoannual

    pseudo[:LongRunInflation].name = "Long Run Inflation"
    pseudo[:LongRunInflation].longname = "Long Run Inflation"
    pseudo[:LongRunInflation].rev_transform = quartertoannual

    pseudo[:MarginalCost].name = "Marginal Cost"
    pseudo[:MarginalCost].longname = "Marginal Cost"

    pseudo[:Wages].name = "Wages"
    pseudo[:Wages].longname = "Wages"

    pseudo[:FlexibleWages].name = "Flexible Wages"
    pseudo[:FlexibleWages].longname = "Wages that would prevail in a flexible-wage economy"

    pseudo[:Hours].name = "Hours"
    pseudo[:Hours].longname = "Hours"

    pseudo[:FlexibleHours].name     = "Flexible Hours"
    pseudo[:FlexibleHours].longname = "Flexible Hours"

    pseudo[:z_t].name     = "z_t"
    pseudo[:z_t].longname = "z_t"

    pseudo[:Expected10YearRateGap].name     = "Expected 10-Year Rate Gap"
    pseudo[:Expected10YearRateGap].longname = "Expected 10-Year Rate Gap"
    pseudo[:Expected10YearRateGap].rev_transform = quartertoannual

    pseudo[:NominalFFR].name     = "Nominal FFR"
    pseudo[:NominalFFR].longname = "Nominal FFR at an annual rate"
    pseudo[:NominalFFR].rev_transform = quartertoannual

    pseudo[:Expected10YearRate].name     = "Expected 10-Year Rate"
    pseudo[:Expected10YearRate].longname = "Expected 10-Year Interest Rate"
    pseudo[:Expected10YearRate].rev_transform = quartertoannual

    pseudo[:Expected10YearNaturalRate].name     = "Expected 10-Year Natural Rate"
    pseudo[:Expected10YearNaturalRate].longname = "Expected 10-Year Natural Rate of Interest"
    pseudo[:Expected10YearNaturalRate].rev_transform = quartertoannual

    pseudo[:ExpectedNominalNaturalRate].name     = "Expected Nominal Natural Rate"
    pseudo[:ExpectedNominalNaturalRate].longname = "Natural Rate + Expected Inflation"
    pseudo[:ExpectedNominalNaturalRate].rev_transform = quartertoannual

    pseudo[:NominalRateGap].name     = "Nominal Rate Gap"
    pseudo[:NominalRateGap].longname = "Nominal FFR - Nominal Natural Rate"
    pseudo[:NominalRateGap].rev_transform = quartertoannual

    pseudo[:LaborProductivityGrowth].name     = "Labor Productivity Growth"
    pseudo[:LaborProductivityGrowth].longname = "Labor Productivity Growth Rate"
    pseudo[:LaborProductivityGrowth].rev_transform = quartertoannual

    # Add to model object
    m.pseudo_observable_mappings = pseudo
end