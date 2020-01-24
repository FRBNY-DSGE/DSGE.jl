function init_pseudo_observable_mappings!(m::SmetsWoutersOrig)

    pseudo_names = [:y_t, :MarginalCost, :Wages, :Hours, :z_t,
                    :NominalFFR, :NominalWageGrowth, :LaborProductivityGrowthNoME,
                    :laborshare_t]

    # Create PseudoObservable objects
    pseudo = OrderedDict{Symbol,PseudoObservable}()
    for k in pseudo_names
        pseudo[k] = PseudoObservable(k)
    end

    pseudo[:y_t].name = "Output Growth"
    pseudo[:y_t].longname = "Output Growth Per Capita"

    pseudo[:MarginalCost].name = "Marginal Cost"
    pseudo[:MarginalCost].longname = "Marginal Cost"

    pseudo[:Wages].name = "Wages"
    pseudo[:Wages].longname = "Wages"

    pseudo[:Hours].name = "Hours"
    pseudo[:Hours].longname = "Hours"

    pseudo[:z_t].name     = "z_t (Technology Growth minus Steady State Growth)"
    pseudo[:z_t].longname = "z_t (Technology Growth minus Steady State Growth)"

    pseudo[:NominalFFR].name     = "Nominal FFR"
    pseudo[:NominalFFR].longname = "Nominal FFR at an annual rate"
    pseudo[:NominalFFR].rev_transform = quartertoannual

    pseudo[:NominalWageGrowth].name = "Nominal Wage Growth"
    pseudo[:NominalWageGrowth].longname = "Nominal Wage Growth"

    pseudo[:LaborProductivityGrowthNoME].name     = "Labor Productivity Growth (No ME)"
    pseudo[:LaborProductivityGrowthNoME].longname = "Labor Productivity Growth (No ME)"

    pseudo[:laborshare_t].name     = "Log Labor Share"
    pseudo[:laborshare_t].longname = "Log Labor Share"

    m.pseudo_observable_mappings = pseudo
end
