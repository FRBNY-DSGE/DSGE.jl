using DSGE, Iterators, Plots

function hair_plot(var::Symbol, df::DataFrame,
                   histories::Vector{MeansBands}, forecasts::Vector{MeansBands};
                   output_file::AbstractString = "",
                   hist_label::AbstractString = "Realized",
                   forecast_label::AbstractString = "Forecasts",
                   forecast_color::Colorant = RGBA(1., 0., 0., 1.))
    # Dates
    dates      = map(quarter_date_to_number, df[:date])
    start_date = ceil(dates[1] / 5) * 5
    end_date   = dates[end]
    date_ticks = start_date:5:end_date

    # Initialize GR backend
    gr()
    p = Plots.plot(xtick = date_ticks)

    # Plot realized (transformed) series
    plot!(p, dates, df[var], label = hist_label, linewidth = 2, linecolor = :black)

    # Plot each forecast
    for (history, forecast) in zip(histories, forecasts)
        dates = vcat([history.means[end, :date]], forecast.means[:date])
        dates = map(quarter_date_to_number, dates)
        series = vcat([history.means[end, var]], forecast.means[var])
        label = history == histories[1] ? forecast_label : ""
        plot!(p, dates, series, label = label, linewidth = 1, linecolor = forecast_color)
    end

    # Save if `output_file` provided
    if !isempty(output_file)
        output_dir = dirname(output_file)
        !isdir(output_dir) && mkdir(output_dir)
        Plots.savefig(output_file)
        println("Saved $output_file")
    end

    return p
end
