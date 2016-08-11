"""
```
transform_data(m::AbstractModel, levels::DataFrame; verbose::Symbol = :low)
```

Transform data loaded in levels and order columns appropriately for the DSGE model. Returns
DataFrame of transformed data.

The DataFrame `levels` is output from `load_data_levels`. The series in levels are
transformed as specified in `m.data_transforms`.
- To prepare for per-capita transformations, population data are filtered using
    `hpfilter`. The series in `levels` to use as the population series is given by the
    `population_mnemonic` setting. If `use_population_forecast` is `true`, a population
    forecast is appended to the recorded population levels before the filtering. Both
    filtered and unfiltered population levels and growth rates are added to the `levels`
    data frame.
- The transformations are applied for each series using the `levels` DataFrame as input.
"""
function transform_data(m::AbstractModel, levels::DataFrame; cond_type::Symbol = :none, verbose::Symbol = :low)

    population_mnemonic = get_setting(m, :population_mnemonic)
    n_obs, _ = size(levels)

    # Step 1: HP filter population forecasts, if they're being used

    # population_recorded: historical population, unfiltered
    # population_all: full unfiltered series (including forecast)
    # dlpopulation_forecast: growth rates of population forecasts pre-filtering

    population_recorded = levels[:,[:date, population_mnemonic]]
    population_all, dlpopulation_forecast, n_population_forecast_obs = if use_population_forecast(m)
        if VERBOSITY[verbose] >= VERBOSITY[:high]
            println("Loading population forecast...")
        end
        
        # load population forecast
        population_forecast_file = inpath(m, "data", "population_forecast_$(data_vintage(m)).csv")
        pop_forecast = readtable(population_forecast_file)

        rename!(pop_forecast, :POPULATION,  population_mnemonic)
        DSGE.na2nan!(pop_forecast)
        DSGE.format_dates!(:date, pop_forecast)

        if cond_type in [:semi, :full]
            pop_forecast = pop_forecast[2:end, :]
        end

        # use our "real" series as current value
        pop_all = vcat(population_recorded, pop_forecast[2:end, :])

        # return values
        pop_all[population_mnemonic],
        difflog(pop_forecast[population_mnemonic]),
        length(pop_forecast[population_mnemonic])
    else
        population_recorded[:,population_mnemonic], [NaN], 1
    end

    # hp filter
    population_all = convert(Array, population_all)
    filtered_population, _ = hpfilter(population_all, 1600)

    # filtered series (levels)
    filtered_population_recorded = filtered_population[1:end-n_population_forecast_obs+1]
    filtered_population_forecast = filtered_population[end-n_population_forecast_obs+1:end]

    # filtered growth rates 
    dlpopulation_recorded          = difflog(population_recorded[population_mnemonic])
    dlfiltered_population_recorded = difflog(filtered_population_recorded)

    levels[:filtered_population]          = filtered_population_recorded
    levels[:filtered_population_growth]   = dlfiltered_population_recorded
    levels[:unfiltered_population_growth] = dlpopulation_recorded

    # Step 2: apply transformations to each series
    transformed = DataFrame()
    transformed[:date] = levels[:date]

    observables = if cond_type == :none
        keys(m.data_transforms)
    elseif cond_type == :semi
        get_setting(m, :cond_semi_names)
    elseif cond_type == :full
        get_setting(m, :cond_full_names)
    end

    for series in observables
        if VERBOSITY[verbose] >= VERBOSITY[:high]
            println("Transforming series $series...")
        end
        f = m.data_transforms[series]
        transformed[series] = f(levels)
    end

    sort!(transformed, cols = :date)

    return transformed
end
