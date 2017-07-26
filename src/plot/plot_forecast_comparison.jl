"""
```
plot_forecast_comparison(var, histold, fcastold, histnew, fcastnew;
    output_file = "", start_date = Nullable{Date}(),
    end_date = Nullable{Date}(), bandpcts::Vector{String} = [\"90.0%\"],
    hist_label = \"History\", old_fcast_label = \"Old Forecast\",
    new_fcast_label = \"New Forecast\", hist_color = :black,
    old_fcast_color = :blue, new_fcast_color = :red, tick_size = 2,
    legend = :best)
```

### Inputs

- `var::Symbol`: e.g. `:obs_gdp`
- `histold::MeansBands`
- `fcastold::MeansBands`
- `histnew::MeansBands`
- `fcastnew::MeansBands`

### Keyword Arguments

- `output_file::String`: if specified, plot will be saved there as a PDF
- `start_date::Nullable{Date}`
- `end_date::Nullable{Date}`
- `bandpcts::Vector{String}`: which bands to plot
- `hist_label::String`
- `old_fcast_label::String`
- `new_fcast_label::String`
- `hist_color::Colorant`
- `old_fcast_color::Colorant`
- `new_fcast_color::Colorant`
- `tick_size::Int`: x-axis (time) tick size in units of years
- `legend`

### Output

- `p::Plot`
"""
function plot_forecast_comparison(var::Symbol,
                                  histold::MeansBands, fcastold::MeansBands,
                                  histnew::MeansBands, fcastnew::MeansBands;
                                  output_file::String = "",
                                  start_date::Nullable{Date} = Nullable{Date}(),
                                  end_date::Nullable{Date} = Nullable{Date}(),
                                  bandpcts::Vector{String} = ["90.0%"],
                                  hist_label::String = "History",
                                  old_hist_label::String="",
                                  old_fcast_label::String = "Old Forecast",
                                  new_fcast_label::String = "New Forecast",
                                  hist_color::Colorant = parse(Colorant, :black),
                                  old_hist_color::Colorant = parse(Colorant, :grey),
                                  old_fcast_color::Colorant = parse(Colorant, :blue),
                                  new_fcast_color::Colorant = parse(Colorant, :red),
                                  tick_size::Int = 2,
                                  legend = :best)

    allold = cat(histold, fcastold)
    allnew = cat(histnew, fcastnew)

    # Convert Dates to datenums (see util.jl)
    for mb in [histold, fcastold, allold, histnew, fcastnew, allnew]
        mb.means[:datenum]      = map(quarter_date_to_number, mb.means[:date])
        mb.bands[var][:datenum] = map(quarter_date_to_number, mb.bands[var][:date])
    end

    # Initialize plot
    p = Plots.plot(legend = legend)

    # Plot new data
    start_date_new, end_date_new = get_date_limits(start_date, end_date, allnew.means[:date])
    start_ind,      end_ind      = get_date_limit_indices(start_date_new, end_date_new, allnew.means[:date])

    n_hist_periods = size(histnew.means, 1)
    hist_inds = start_ind:min(end_ind, n_hist_periods)

    plot!(p, histnew.means[hist_inds, :datenum], histnew.means[hist_inds, var],
          linewidth = 2, linecolor = hist_color, label = hist_label)

    # Plot new forecast
    n_hist_periods  = size(histnew.means, 1)
    n_all_periods   = size(allnew.means, 1)
    fcast_inds      = max(start_ind, n_hist_periods):end_ind

    plot!(p, allnew.means[fcast_inds, :datenum], allnew.means[fcast_inds, var],
          linewidth = 2, linecolor = new_fcast_color, label = new_fcast_label)

    if fcastnew.metadata[:para] in [:full, :subset]
        for pct in bandpcts
            plot!(p, allnew.bands[var][fcast_inds, :datenum], allnew.bands[var][fcast_inds, Symbol(pct, " UB")],
                  linewidth = 2, linecolor = new_fcast_color, label = "")
            plot!(p, allnew.bands[var][fcast_inds, :datenum], allnew.bands[var][fcast_inds, Symbol(pct, " LB")],
                  linewidth = 2, linecolor = new_fcast_color, label = "")
        end
    end

    # Plot old data
    # Conditional on whether the historic path differs (i.e. for pseudo-observables)
    start_date_old, end_date_old = get_date_limits(start_date, end_date, allold.means[:date])
    start_ind,      end_ind      = get_date_limit_indices(start_date_old, end_date_old, allold.means[:date])

    if !isempty(old_hist_label)
        n_hist_periods = size(histold.means, 1)
        hist_inds = start_ind:min(end_ind, n_hist_periods)

        plot!(p, histold.means[hist_inds, :datenum], histold.means[hist_inds, var],
              linewidth = 2, linecolor = old_hist_color, label = old_hist_label)
    end

    # Plot old forecast
    n_hist_periods  = size(histold.means, 1)
    n_all_periods   = size(allold.means, 1)
    fcast_inds      = max(start_ind, n_hist_periods):end_ind

    plot!(p, allold.means[fcast_inds, :datenum], allold.means[fcast_inds, var],
          linewidth = 2, linecolor = old_fcast_color, linestyle = :dash, label = old_fcast_label)

    if fcastold.metadata[:para] in [:full, :subset]
        for pct in bandpcts
            plot!(p, allold.bands[var][fcast_inds, :datenum], allold.bands[var][fcast_inds, Symbol(pct, " UB")],
                  linewidth = 2, linecolor = old_fcast_color, linestyle = :dash, label = "")
            plot!(p, allold.bands[var][fcast_inds, :datenum], allold.bands[var][fcast_inds, Symbol(pct, " LB")],
                  linewidth = 2, linecolor = old_fcast_color, linestyle = :dash, label = "")
        end
    end

    # Set date ticks
    date_ticks!(p, min(start_date_old, start_date_new), max(end_date_old, end_date_new), tick_size)

    # Save if output_file provided
    save_plot(p, output_file)

    # Delete datenum fields from input arguments (causes problems in cat)
    for mb in [histold, fcastold, histnew, fcastnew]
        delete!(mb.means,      :datenum)
        delete!(mb.bands[var], :datenum)
    end

    return p
end
