"""
`transform_data(m::AbstractModel, levels::DataFrame)`

Transform data loaded in levels and order columns appropriately for the DSGE model.
"""
function transform_data(m::AbstractModel, levels::DataFrame)

    n_obs, _ = size(levels)

    # HP filter population forecasts, if they're being used

    # popreal: historical population, unfiltered
    # population_all: filtered full series (including Macro Advisers forecast)
    # dlpopforecast: growth rates of population forecasts pre-filtering

    popreal = levels[:,[:date, :CNP16OV]]
    population_all, dlpopforecast, n_popforecast_obs = if use_population_forecast(m)

        # load population forecast

        #TODO: refactor "popforecast" etc. -> "population_forecast" everywhere
        popforecast_file = inpath(m, "data", "popforecast_$(data_vintage(m)).txt")
        #TODO: please use CSV :)
        pop_forecast = readtable(popforecast_file, separator = '\t')
        #TODO: confused why CNP16OV is hardcoded here
        rename!(pop_forecast, :POPULATION,  :CNP16OV)
        DSGE.na2nan!(pop_forecast)
        DSGE.format_dates!(:date, pop_forecast)

        # use our "real" series as current value
        pop_all = [popreal; pop_forecast[2:end,:]]

        pop_all[:CNP16OV], difflog(pop_forecast[:CNP16OV]), length(pop_forecast[:CNP16OV])
    else
        popreal, _, 0
    end

    # hp filter
    population_all = convert(Array, population_all)
    popfor, _ = hpfilter(population_all, 1600)

    # filtered series (levels)
    population = popfor[1:end-n_popforecast_obs+1]     # filtered recorded population series
    MA_pop     = popfor[end-n_popforecast_obs+1:end]   # filtered forecast population series

    # growth rates - TODO: clean up, leaving all names now for easier comparison to MATLAB
    dlpopreal = difflog(popreal[:CNP16OV])
    dlpopall = difflog(population)

    levels[:filtered_population] = population
    levels[:filtered_population_growth] = dlpopall
    levels[:unfiltered_population_growth] = dlpopreal

    # Now apply transformations to each series
    transformed = DataFrame()
    transformed[:date] = levels[:date]

    #TODO can we remove enumerate?
    for series in keys(m.data_transforms)
        f = m.data_transforms[series]
        transformed[series] = f(levels)
    end

    sort!(transformed, cols = :date)
end
