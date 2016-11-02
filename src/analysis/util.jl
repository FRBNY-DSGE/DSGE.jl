function resize_population_forecast(population_forecast::DataFrame, nperiods::Int;
                                           population_mnemonic::Symbol = Symbol())

    # number of periods to extend population forecast
    n_filler_periods = nperiods - size(population_forecast,1)

    # Extract population mnemonic from Nullable object (if not null). If null,
    # take a guess or throw an error if you can't tell.
    mnemonic = if population_mnemonic == Symbol()
        if size(population_forecast,2) > 2
            error("Please indicate which column contains population
                forecasts using the population_mnemonic keyword argument")
        end

        setdiff(names(population_forecast), [:date])[1]
    else
        population_mnemonic
    end

    last_provided = population_forecast[end,:date]

    # create date range. There are on average 91.25 days in a quarter.
    dr = last_provided:(last_provided+Dates.Day(93 * n_filler_periods))

    islastdayofquarter = x->Dates.lastdayofquarter(x) == x
    dates = recur(dr) do x
        islastdayofquarter(x)
    end

    # first element of dates is the last quarter of population forecasts supplied, so we
    # cut it off here
    dates = dates[2:n_filler_periods+1]

    extra = DataFrame()
    extra[:date] = dates

    # resize population forecast by adding n_filler_periods to the given forecast.
    resized = if n_filler_periods > 0

        extra_stuff = fill(population_forecast[end,mnemonic], n_filler_periods)
        extra[mnemonic] = extra_stuff

        [population_forecast; extra]
    elseif n_filler_periods < 0
        population_forecast[1:nperiods,:]
    else
        population_forecast
    end

    resized
end