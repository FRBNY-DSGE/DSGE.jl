function init_pseudo_observable_mappings!(m::Model1010)

    pseudo_names = if subspec(m) in ["ss2", "ss4"]
	    [:y_t, :y_f_t, :OutputGap,
	     :π_t, :LongRunInflation,
	     :Wages, :FlexibleWages, :z_t, :Hours, :FlexibleHours,
	     :RealNaturalRate, :ExAnteRealRate,
	     :NominalFFR, :ExpectedAvgNominalNaturalRate, :NominalRateGap,
	     :ExpectedAvg10YearRealRate,    :ExpectedAvg10YearRealNaturalRate,
	     :ExpectedAvg10YearNominalRate, :ExpectedAvg10YearNominalNaturalRate,:ExpectedAvg10YearRateGap,
	     :ExpectedAvg5YearNominalRate,  :ExpectedAvg5YearNominalNaturalRate,
	     :ExpectedAvg5YearRealRate,     :ExpectedAvg5YearRealNaturalRate, :ExpectedAvg5YearRateGap,
	     :ExpectedAvg20YearNominalRate, :ExpectedAvg20YearNominalNaturalRate,
	     :ExpectedAvg20YearRealRate,    :ExpectedAvg20YearRealNaturalRate,   :ExpectedAvg20YearRateGap,
	     :Forward5YearRateGap,
	     :Forward10YearRateGap,
	     :Forward20YearRealRate, :Forward20YearRealNaturalRate,
	     :Forward30YearRealRate, :Forward30YearRealNaturalRate]
    else
        [:y_t, :y_f_t, :OutputGap,
         :Hours,
         :RealNaturalRate, :ExAnteRealRate,
         :NominalNaturalRate,
         :NominalFFR,
         :ExpectedAvg5YearRealRate, :ExpectedAvg5YearRealNaturalRate,
         :ExpectedAvg10YearRealRate,:ExpectedAvg10YearRealNaturalRate,
         :ExpectedAvg20YearRealRate, :ExpectedAvg20YearRealNaturalRate,
         :Forward5YearRealRate,  :Forward5YearRealNaturalRate,
         :Forward10YearRealRate, :Forward10YearRealNaturalRate,
         :Forward20YearRealRate, :Forward20YearRealNaturalRate,
         :Forward30YearRealRate, :Forward30YearRealNaturalRate]
    end

    if subspec(m) in ["ss13", "ss18"]
        push!(pseudo_names, :LiquidityConvenienceYield, :SafetyConvenienceYield,
              :Forward20YearLiquidityConvenienceYield, :Forward20YearSafetyConvenienceYield,
              :SDF, :Forward20YearSDF)
    end

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

	pseudo[:OutputGap].name = "Output Gap"
	pseudo[:OutputGap].longname = "Output Gap"

    if subspec(m) in ["ss2", "ss4"]
	    pseudo[:π_t].name = "Inflation"
	    pseudo[:π_t].longname = "Inflation"
	    pseudo[:π_t].rev_transform = quartertoannual

	    pseudo[:LongRunInflation].name = "Long Run Inflation"
	    pseudo[:LongRunInflation].longname = "Long Run Inflation"
	    pseudo[:LongRunInflation].rev_transform = quartertoannual

	    pseudo[:Wages].name = "Wages"
	    pseudo[:Wages].longname = "Wages"

	    pseudo[:FlexibleWages].name = "Flexible Wages"
	    pseudo[:FlexibleWages].longname = "Wages that would prevail in a flexible-wage economy"

	    pseudo[:z_t].name     = "z_t"
	    pseudo[:z_t].longname = "z_t"
    end

	pseudo[:Hours].name = "Hours"
	pseudo[:Hours].longname = "Hours"

    if subspec(m) in ["ss2", "ss4"]
	    pseudo[:FlexibleHours].name     = "Flexible Hours"
	    pseudo[:FlexibleHours].longname = "Flexible Hours"
    end

	pseudo[:RealNaturalRate].name = "Real Natural Rate"
	pseudo[:RealNaturalRate].longname = "The real interest rate that would prevail in a flexible-price economy."
	pseudo[:RealNaturalRate].rev_transform = quartertoannual

	pseudo[:ExAnteRealRate].name = "Ex Ante Real Rate"
	pseudo[:ExAnteRealRate].longname = "Ex Ante Real Rate"
	pseudo[:ExAnteRealRate].rev_transform = quartertoannual

	pseudo[:NominalFFR].name     = "Nominal FFR"
	pseudo[:NominalFFR].longname = "Nominal FFR at an annual rate"
	pseudo[:NominalFFR].rev_transform = quartertoannual

    if subspec(m) in ["ss2", "ss4"]
        pseudo[:ExpectedAvgNominalNaturalRate].name     = "Expected Average Nominal Natural Rate"
        pseudo[:ExpectedAvgNominalNaturalRate].longname = "Natural Rate + Expected Inflation"
        pseudo[:ExpectedAvgNominalNaturalRate].rev_transform = quartertoannual

        pseudo[:NominalRateGap].name     = "Nominal Rate Gap"
        pseudo[:NominalRateGap].longname = "Nominal FFR - Nominal Natural Rate"
        pseudo[:NominalRateGap].rev_transform = quartertoannual
    end

    pseudo[:ExpectedAvg10YearRealRate].name     = "Expected Average 10-Year Real Rate"
    pseudo[:ExpectedAvg10YearRealRate].longname = "Expected Average 10-Year Real Interest Rate"
    pseudo[:ExpectedAvg10YearRealRate].rev_transform = quartertoannual

    pseudo[:ExpectedAvg10YearRealNaturalRate].name     = "Expected Average 10-Year Real Natural Rate"
    pseudo[:ExpectedAvg10YearRealNaturalRate].longname = "Expected Average 10-Year Real Natural Rate of Interest"
    pseudo[:ExpectedAvg10YearRealNaturalRate].rev_transform = quartertoannual

    if subspec(m) in ["ss2", "ss4"]
        pseudo[:ExpectedAvg10YearNominalRate].name     = "Expected Average 10-Year Nominal Rate"
        pseudo[:ExpectedAvg10YearNominalRate].longname = "Expected Average 10-Year Nominal Interest Rate"
        pseudo[:ExpectedAvg10YearNominalRate].rev_transform = quartertoannual

        pseudo[:ExpectedAvg10YearNominalNaturalRate].name     = "Expected Average 10-Year Nominal Natural Rate"
        pseudo[:ExpectedAvg10YearNominalNaturalRate].longname = "Expected Average 10-Year Nominal Natural Rate of Interest"
        pseudo[:ExpectedAvg10YearNominalNaturalRate].rev_transform = quartertoannual

        pseudo[:ExpectedAvg10YearRateGap].name     = "Expected Average 10-Year Rate Gap"
        pseudo[:ExpectedAvg10YearRateGap].longname = "Expected Average 10-Year Rate Gap"
        pseudo[:ExpectedAvg10YearRateGap].rev_transform = quartertoannual
    end

    pseudo[:ExpectedAvg5YearRealRate].name     = "Expected Average 5-Year Real Rate"
    pseudo[:ExpectedAvg5YearRealRate].longname = "Expected Average 5-Year Real Interest Rate"
    pseudo[:ExpectedAvg5YearRealRate].rev_transform = quartertoannual

    pseudo[:ExpectedAvg5YearRealNaturalRate].name     = "Expected Average 5-Year Real Natural Rate"
    pseudo[:ExpectedAvg5YearRealNaturalRate].longname = "Expected Average 5-Year Real Natural Rate of Interest"
    pseudo[:ExpectedAvg5YearRealNaturalRate].rev_transform = quartertoannual

    if subspec(m) in ["ss2", "ss4"]
        pseudo[:ExpectedAvg5YearNominalRate].name     = "Expected Average 5-Year Rate"
        pseudo[:ExpectedAvg5YearNominalRate].longname = "Expected Average 5-Year Interest Rate"
        pseudo[:ExpectedAvg5YearNominalRate].rev_transform = quartertoannual

        pseudo[:ExpectedAvg5YearNominalNaturalRate].name     = "Expected Average 5-Year Natural Rate"
        pseudo[:ExpectedAvg5YearNominalNaturalRate].longname = "Expected Average 5-Year Natural Rate of Interest"
        pseudo[:ExpectedAvg5YearNominalNaturalRate].rev_transform = quartertoannual

        pseudo[:ExpectedAvg5YearRateGap].name     = "Expected Average 5-Year Rate Gap"
        pseudo[:ExpectedAvg5YearRateGap].longname = "Expected Average 5-Year Rate Gap"
        pseudo[:ExpectedAvg5YearRateGap].rev_transform = quartertoannual
    end

    pseudo[:ExpectedAvg20YearRealRate].name     = "Expected Average 20-Year Real Rate"
    pseudo[:ExpectedAvg20YearRealRate].longname = "Expected Average 20-Year Real Interest Rate"
    pseudo[:ExpectedAvg20YearRealRate].rev_transform = quartertoannual

    pseudo[:ExpectedAvg20YearRealNaturalRate].name     = "Expected Average 20-Year Real Natural Rate"
    pseudo[:ExpectedAvg20YearRealNaturalRate].longname = "Expected Average 20-Year Real Natural Rate of Interest"
    pseudo[:ExpectedAvg20YearRealNaturalRate].rev_transform = quartertoannual

    if subspec(m) in ["ss2", "ss4"]
        pseudo[:ExpectedAvg20YearNominalRate].name     = "Expected Average 20-Year Nominal Rate"
        pseudo[:ExpectedAvg20YearNominalRate].longname = "Expected Average 20-Year Nominal Interest Rate"
        pseudo[:ExpectedAvg20YearNominalRate].rev_transform = quartertoannual

        pseudo[:ExpectedAvg20YearNominalNaturalRate].name     = "Expected Average 20-Year Nominal Natural Rate"
        pseudo[:ExpectedAvg20YearNominalNaturalRate].longname = "Expected Average 20-Year Nominal Natural Rate of Interest"
        pseudo[:ExpectedAvg20YearNominalNaturalRate].rev_transform = quartertoannual

        pseudo[:ExpectedAvg20YearRateGap].name     = "Expected Average 20-Year Rate Gap"
        pseudo[:ExpectedAvg20YearRateGap].longname = "Expected Average 20-Year Rate Gap"
        pseudo[:ExpectedAvg20YearRateGap].rev_transform = quartertoannual

        pseudo[:Forward5YearRateGap].name     = "5-Year Forward Rate Gap"
        pseudo[:Forward5YearRateGap].longname = "5-Year Forward Rate Gap"
        pseudo[:Forward5YearRateGap].rev_transform = quartertoannual

        pseudo[:Forward10YearRateGap].name     = "10-Year Forward Rate Gap"
        pseudo[:Forward10YearRateGap].longname = "10-Year Forward Rate Gap"
        pseudo[:Forward10YearRateGap].rev_transform = quartertoannual
    end

    pseudo[:Forward20YearRealRate].name     = "Forward 20-Year Real Rate"
    pseudo[:Forward20YearRealRate].longname = "Forward 20-Year Real Rate of Interest (not the average, computed by projecting TTT foreward 80 periods)"
    pseudo[:Forward20YearRealRate].rev_transform = quartertoannual

    pseudo[:Forward20YearRealNaturalRate].name     = "Forward 20-Year Real Natural Rate"
    pseudo[:Forward20YearRealNaturalRate].longname = "Forward 20-Year Real Natural Rate of Interest (not the average, computed by projecting TTT foreward 80 periods)"
    pseudo[:Forward20YearRealNaturalRate].rev_transform = quartertoannual

    pseudo[:Forward30YearRealRate].name     = "Forward 30-Year Real Rate"
    pseudo[:Forward30YearRealRate].longname = "Forward 30-Year Real Rate of Interest (not the average, computed by projecting TTT foreward 80 periods)"
    pseudo[:Forward30YearRealRate].rev_transform = quartertoannual

    pseudo[:Forward30YearRealNaturalRate].name     = "Forward 30-Year Real Natural Rate"
    pseudo[:Forward30YearRealNaturalRate].longname = "Forward 30-Year Real Natural Rate of Interest (not the average, computed by projecting TTT foreward 80 periods)"
    pseudo[:Forward30YearRealNaturalRate].rev_transform = quartertoannual

    if !(subspec(m) in ["ss2", "ss4"])
        pseudo[:NominalNaturalRate].name     = "Nominal Natural Rate"
        pseudo[:NominalNaturalRate].longname = "Natural Rate + Expected Inflation"
        pseudo[:NominalNaturalRate].rev_transform = quartertoannual

        pseudo[:Forward5YearRealRate].name     = "Forward 5-Year Real Rate"
        pseudo[:Forward5YearRealRate].longname = "Forward 5-Year Real Rate of Interest (not the average, computed by projecting TTT foreward 80 periods)"
        pseudo[:Forward5YearRealRate].rev_transform = quartertoannual

        pseudo[:Forward5YearRealNaturalRate].name     = "Forward 5-Year Real Natural Rate"
        pseudo[:Forward5YearRealNaturalRate].longname = "Forward 5-Year Real Natural Rate of Interest (not the average, computed by projecting TTT foreward 80 periods)"
        pseudo[:Forward5YearRealNaturalRate].rev_transform = quartertoannual

        pseudo[:Forward10YearRealRate].name     = "Forward 10-Year Real Rate"
        pseudo[:Forward10YearRealRate].longname = "Forward 10-Year Real Rate of Interest (not the average, computed by projecting TTT foreward 80 periods)"
        pseudo[:Forward10YearRealRate].rev_transform = quartertoannual

        pseudo[:Forward10YearRealNaturalRate].name     = "Forward 10-Year Real Natural Rate"
        pseudo[:Forward10YearRealNaturalRate].longname = "Forward 10-Year Real Natural Rate of Interest (not the average, computed by projecting TTT foreward 80 periods)"
        pseudo[:Forward10YearRealNaturalRate].rev_transform = quartertoannual
    end

    if subspec(m) in ["ss13", "ss18"]
        pseudo[:LiquidityConvenienceYield].name = "Liquidity convenience yield"
        pseudo[:LiquidityConvenienceYield].longname = "Liquidity convenience yield"
        pseudo[:LiquidityConvenienceYield].rev_transform = quartertoannual

        pseudo[:SafetyConvenienceYield].name = "Safety convenience yield"
        pseudo[:SafetyConvenienceYield].longname = "Safety convenience yield"
        pseudo[:SafetyConvenienceYield].rev_transform = quartertoannual

        pseudo[:SDF].name = "Stochastic discount factor"
        pseudo[:SDF].longname = "Stochastic discount factor"
        pseudo[:SDF].rev_transform = quartertoannual

        pseudo[:Forward20YearLiquidityConvenienceYield].name = "20-year forward liquidity convenience yield"
        pseudo[:Forward20YearLiquidityConvenienceYield].longname = "20-year forward liquidity convenience yield"
        pseudo[:Forward20YearLiquidityConvenienceYield].rev_transform = quartertoannual

        pseudo[:Forward20YearSafetyConvenienceYield].name = "20-year forward safety convenience yield"
        pseudo[:Forward20YearSafetyConvenienceYield].longname = "20-year forward safety convenience yield"
        pseudo[:Forward20YearSafetyConvenienceYield].rev_transform = quartertoannual

        pseudo[:Forward20YearSDF].name = "20-year forward stochastic discount factor"
        pseudo[:Forward20YearSDF].longname = "20-year forward stochastic discount factor"
        pseudo[:Forward20YearSDF].rev_transform = quartertoannual
    end

    # Add to model object
    m.pseudo_observable_mappings = pseudo
end
