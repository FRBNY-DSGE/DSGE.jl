"""
```
plot_history_and_forecast(var, history, forecast; output_file = "",
    start_date = Nullable{Date}(), end_date = Nullable{Date}(),
    hist_label = \"History\", forecast_label = \"Forecast\",
    hist_color = :black, forecast_mean_color = :red,
    forecast_band_color = RGBA(0, 0, 1, 0.1), tick_size = 5, legend = :best,
    plot_handle = plot())
```

Plot `var` from `history` and `forecast`. If these correspond to a
full-distribution forecast, the forecast will be a fan chart.

### Inputs

- `var::Symbol`: e.g. `:obs_gdp`
- `history::MeansBands` or `Vector{MeansBands}`
- `forecast::MeansBands` or `Vector{MeansBands}`

### Keyword Arguments

- `output_file::String`: if specified, plot will be saved there as a PDF
- `start_date::Nullable{Date}`
- `end_date::Nullable{Date}`
- `hist_label::String`
- `forecast_label::String`
- `hist_color::Colorant`
- `forecast_mean_color::Colorant`
- `forecast_band_color::Colorant`
- `tick_size::Int`: x-axis (time) tick size in units of years
- `legend`
- `plot_handle::Plots.Plot`: a plot handle to add `history` and `forecast` to

### Output

- `p::Plot`
"""
function plot_history_and_forecast(var::Symbol, history::MeansBands, forecast::MeansBands;
                                   output_file::String = "",
                                   start_date::Nullable{Date} = Nullable{Date}(),
                                   end_date::Nullable{Date} = Nullable{Date}(),
                                   hist_label::String = "History",
                                   forecast_label::String = "Forecast",
                                   hist_color::Colorant = RGBA(0., 0., 0., 1.),
                                   forecast_mean_color::Colorant = RGBA(1., 0., 0., 1.),
                                   forecast_band_color::Colorant = RGBA(0., 0., 1., 0.1),
                                   tick_size::Int = 5,
                                   legend = :best,
                                   plot_handle::Plots.Plot = plot())
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

    # Initialize plot
    p = if isempty(plot_handle.series_list)
        Plots.plot(legend = legend)
    else
        plot_handle
    end

    # Plot bands
    if combined.metadata[:para] in [:full, :subset]
        band_inds = get_bands_indices(var, history, forecast, hist_inds, fcast_inds)
        band_pcts = DSGE.which_density_bands(combined; uniquify = true)
        for pct in band_pcts
            plot!(p, datenums[band_inds], combined.bands[var][band_inds, Symbol("$pct UB")],
                  fillto = combined.bands[var][band_inds, Symbol("$pct LB")],
                  label = "", color = forecast_band_color, α = 0.10)
        end
    end

    # Plot mean
    plot!(p, datenums[hist_inds],  combined.means[hist_inds,  var], label = hist_label,
          linewidth = 2, linecolor = hist_color)
    plot!(p, datenums[fcast_inds], combined.means[fcast_inds, var], label = forecast_label,
          linewidth = 2, linecolor = forecast_mean_color)

    # Set date ticks
    date_ticks!(p, start_date, end_date, tick_size)

    # Save if output_file provided
    save_plot(p, output_file)

    return p
end


function plot_history_and_forecast(var::Symbol, histories::Vector{MeansBands}, forecasts::Vector{MeansBands};
                                   start_date::Nullable{Date} = Nullable{Date}(),
                                   end_date::Nullable{Date} = Nullable{Date}(),
                                   output_file::String = "",
                                   hist_label::String = "History",
                                   forecast_label::String = "Forecast",
                                   hist_color::Colorant = RGBA(0., 0., 0., 1.),
                                   forecast_mean_color::Colorant = RGBA(1., 0., 0., 1.),
                                   forecast_band_color::Colorant = RGBA(0., 0., 1., 0.1),
                                   tick_size::Int = 5,
                                   legend = :best,
                                   plot_handle::Plots.Plot = plot())

    @assert length(histories) == length(forecasts) "histories and forecasts must be same length"

    for (history, forecast) in zip(histories, forecasts)
       plot_handle =  plot_history_and_forecast(var, history, forecast;
                            start_date = start_date, end_date = end_date,
                            output_file = output_file,
                            hist_label = hist_label, forecast_label = forecast_label,
                            hist_color = hist_color, forecast_mean_color = forecast_mean_color,
                            forecast_band_color = forecast_band_color,
                            tick_size = tick_size, legend = legend, plot_handle = plot_handle)
    end

    return plot_handle
end