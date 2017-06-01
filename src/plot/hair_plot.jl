using DSGE, Iterators, Plots

function hair_plot(var::Symbol, df::DataFrame,
                   histories::Vector{MeansBands}, forecasts::Vector{MeansBands};
                   output_file::String = "",
                   hist_label::String = "Realized",
                   forecast_label::String = "Forecasts",
                   forecast_palette::Symbol = Symbol(),
                   forecast_color::Colorant = RGBA(1., 0., 0., 1.),
                   legend = :best)

    initial_values = map(history -> history.means[end, var], histories)
    hair_plot(var, df, initial_values, forecasts,
              output_file = output_file, hist_label = hist_label, forecast_label = forecast_label,
              forecast_palette = forecast_palette, forecast_color = forecast_color, legend = legend)
end

function hair_plot(var::Symbol, df::DataFrame,
                   initial_values::Vector{Float64}, forecasts::Vector{MeansBands};
                   output_file::String = "",
                   hist_label::String = "Realized",
                   forecast_label::String = "Forecasts",
                   forecast_palette::Symbol = Symbol(),
                   forecast_color::Colorant = RGBA(1., 0., 0., 1.),
                   legend = :best)
    # Dates
    datenums   = map(quarter_date_to_number, df[:date])
    date_ticks = get_date_ticks(df[:date])

    # Initialize GR backend
    gr()
    p = Plots.plot(xtick = date_ticks, legend = legend)

    # Plot realized (transformed) series
    plot!(p, datenums, df[var], label = hist_label, linewidth = 2, linecolor = :black)

    # Plot each forecast
    for (initial_value, forecast) in zip(initial_values, forecasts)
        date_0 = DSGE.iterate_quarters(forecast.means[1, :date], -1)
        dates = vcat([date_0], forecast.means[:date])
        datenums = map(quarter_date_to_number, dates)

        series = vcat([initial_value], forecast.means[var])

        label = forecast == forecasts[1] ? forecast_label : ""
        if forecast_palette == Symbol()
            plot!(p, datenums, series, label = label, linewidth = 1, linecolor = forecast_color)
        else
            plot!(p, datenums, series, label = label, linewidth = 1, palette = forecast_palette)
        end
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