using DSGE, GR, Plots

function plot_history_and_forecast(var::Symbol, history::MeansBands, forecast::MeansBands;
                                   # start_date::Nullable{Date} = Nullable{Date}(),
                                   # end_date::Nullable{Date} = Nullable{Date}(),
                                   output_file::AbstractString = "",
                                   hist_label::AbstractString = "History",
                                   forecast_label::AbstractString = "Forecast",
                                   hist_color::Colorant = RGBA(0., 0., 0., 1.),
                                   forecast_mean_color::Colorant = RGBA(1., 0., 0., 1.),
                                   forecast_band_color::Colorant = RGBA(0., 0., 1., 0.1))
    # Concatenate MeansBands
    combined = cat(history, forecast)
    hist_inds  = 1:size(history.means, 1)
    fcast_inds = size(history.means, 1):size(combined.means, 1)

    # Dates
    dates      = map(quarter_date_to_number, combined.means[:date])
    start_date = ceil(dates[1] / 5) * 5
    end_date   = dates[end]
    date_ticks = start_date:5:end_date

    # Initialize GR backend
    gr()
    p = Plots.plot(xtick = date_ticks)

    # Plot bands
    band_percents = DSGE.which_density_bands(combined; uniquify = true)
    for pct in band_percents
        plot!(p, dates, combined.bands[var][Symbol("$pct UB")],
              fillto = combined.bands[var][Symbol("$pct LB")],
              label = "", color = forecast_band_color)
              # label = "", color = :blue, α = 0.10)
    end

    # Plot mean
    plot!(p, dates[hist_inds],  combined.means[hist_inds,  var], label = hist_label,
          linewidth = 2, linecolor = hist_color)
    plot!(p, dates[fcast_inds], combined.means[fcast_inds, var], label = forecast_label,
          linewidth = 2, linecolor = forecast_mean_color)

    # Save if `output_file` provided
    if !isempty(output_file)
        output_dir = dirname(output_file)
        !isdir(output_dir) && mkdir(output_dir)
        Plots.savefig(output_file)
        println("Saved $output_file")
    end

    return p
end