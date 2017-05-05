using DSGE, Plots

function plot_history_and_forecast(var::Symbol, history::MeansBands, forecast::MeansBands;
                                   start_date::Nullable{Date} = Nullable{Date}(),
                                   end_date::Nullable{Date} = Nullable{Date}(),
                                   output_file::AbstractString = "",
                                   hist_label::AbstractString = "History",
                                   forecast_label::AbstractString = "Forecast",
                                   hist_color::Colorant = RGBA(0., 0., 0., 1.),
                                   forecast_mean_color::Colorant = RGBA(1., 0., 0., 1.),
                                   forecast_band_color::Colorant = RGBA(0., 0., 1., 0.1),
                                   tick_size::Int = 5,
                                   legend = :best)
    # Concatenate MeansBands
    combined = cat(history, forecast)

    # Dates
    start_date, end_date = get_date_limits(start_date, end_date, combined.means[:date])
    start_ind,  end_ind  = get_date_limit_indices(start_date, end_date, combined.means[:date])
    datenums             = map(quarter_date_to_number, combined.means[:date])

    # Indices
    n_hist_periods  = size(history.means,  1)
    n_fcast_periods = size(forecast.means, 1)
    n_all_periods   = n_hist_periods + n_fcast_periods

    hist_inds  = start_ind:min(end_ind, n_hist_periods)
    fcast_inds = max(start_ind, n_hist_periods):end_ind
    all_inds   = start_ind:end_ind

    # Initialize GR backend
    gr()
    p = Plots.plot(legend = legend)

    # Plot bands
    band_percents = DSGE.which_density_bands(combined; uniquify = true)
    for pct in band_percents
        plot!(p, datenums[all_inds], combined.bands[var][all_inds, Symbol("$pct UB")],
              fillto = combined.bands[var][all_inds, Symbol("$pct LB")],
              label = "", color = forecast_band_color, α = 0.10)
    end

    # Plot mean
    plot!(p, datenums[hist_inds],  combined.means[hist_inds,  var], label = hist_label,
          linewidth = 2, linecolor = hist_color)
    plot!(p, datenums[fcast_inds], combined.means[fcast_inds, var], label = forecast_label,
          linewidth = 2, linecolor = forecast_mean_color)

    # Set x-axis limits and ticks
    # xlims attribute only sets finite values, e.g. (-Inf, 2) sets only the right limit
    t0 = isnull(start_date) ? -Inf : quarter_date_to_number(get(start_date))
    t1 = isnull(end_date)   ?  Inf : quarter_date_to_number(get(end_date))
    date_ticks = get_date_ticks(get(start_date), get(end_date), tick_size = tick_size)
    xaxis!(p, xlims = (t0, t1), xtick = date_ticks)

    # Save if output_file provided
    if !isempty(output_file)
        output_dir = dirname(output_file)
        !isdir(output_dir) && mkdir(output_dir)
        Plots.savefig(output_file)
        println("Saved $output_file")
    end

    return p
end