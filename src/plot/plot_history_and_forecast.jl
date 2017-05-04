using DSGE, Plots

function plot_history_and_forecast(var::Symbol, history::MeansBands, forecast::MeansBands;
                                   start_date::Nullable{Date} = Nullable{Date}(),
                                   output_file::AbstractString = "",
                                   hist_label::AbstractString = "History",
                                   forecast_label::AbstractString = "Forecast",
                                   hist_color::Colorant = RGBA(0., 0., 0., 1.),
                                   forecast_mean_color::Colorant = RGBA(1., 0., 0., 1.),
                                   forecast_band_color::Colorant = RGBA(0., 0., 1., 0.1))
    # Concatenate MeansBands
    combined = cat(history, forecast)

    # Dates
    datenums   = map(quarter_date_to_number, df[:date])
    start_ind  = if !isnull(start_date)
        @assert history.means[1, :date] <= get(start_date) <= history.means[end, :date]
        findfirst(combined.means[:date], get(start_date))
    else
        1
    end
    date_ticks = get_date_ticks(combined.means[start_ind:end, :date])

    hist_inds  = start_ind:size(history.means, 1)
    fcast_inds = size(history.means, 1):size(combined.means, 1)
    all_inds   = start_ind:size(combined.means, 1)

   # Initialize GR backend
    gr()
    p = Plots.plot(xtick = date_ticks)

    # Plot bands
    band_percents = DSGE.which_density_bands(combined; uniquify = true)
    for pct in band_percents
        plot!(p, datenums[all_inds], combined.bands[var][all_inds, Symbol("$pct UB")],
              fillto = combined.bands[var][all_inds, Symbol("$pct LB")],
              label = "", color = forecast_band_color)
              # label = "", color = :blue, α = 0.10)
    end

    # Plot mean
    plot!(p, datenums[hist_inds],  combined.means[hist_inds,  var], label = hist_label,
          linewidth = 2, linecolor = hist_color)
    plot!(p, datenums[fcast_inds], combined.means[fcast_inds, var], label = forecast_label,
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