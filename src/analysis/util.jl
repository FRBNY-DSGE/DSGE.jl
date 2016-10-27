function resize_population_forecast(population_forecast, nperiods)

    n_filler_periods = nperiods - length(population_forecast)

    resized = if n_filler_periods > 0
        [population_forecast; fill(population_forecast[end], n_filler_periods)]
    elseif n_filler_periods < 0
        population_forecast[1:-n_filler_periods]
    else
        population_forecast
    end

    resized
end