"""
```
plot_history_and_forecast(m, var, class, input_type, cond_type;
    forecast_string = "", bdd_and_unbdd = false, output_file = "",
    title = "", kwargs...)

plot_history_and_forecast(m, vars, class, input_type, cond_type;
    forecast_string = "", bdd_and_unbdd = false, output_files = [],
    titles = [], kwargs...)

plot_history_and_forecast(var, history, forecast; output_file = "",
    title = "", start_date = Nullable{Date}(), end_date = Nullable{Date}(),
    hist_label = \"History\", forecast_label = \"Forecast\",
    hist_color = :black, forecast_color = :red, linestyle = :solid,
    bands_color = RGBA(0, 0, 1, 0.1), bands_pcts = [],
    bands_style = :fan, tick_size = 5, legend = :best,
    plot_handle = plot())
```

Plot `var` from `history` and `forecast`, possibly read in using `read_mb`
(depending on the method). If these correspond to a full-distribution forecast,
you can specify the `bands_style` and `bands_pcts`.

### Inputs

- `var::Symbol` or `vars::Vector{Symbol}`: variable(s) to be plotted,
  e.g. `:obs_gdp` or `[:obs_gdp, :obs_nominalrate]`

**Methods 1 and 2 only:**

- `m::AbstractModel`
- `class::Symbol`

**Method 3 only:**

- `history::MeansBands` or `Vector{MeansBands}`
- `forecast::MeansBands` or `Vector{MeansBands}`

### Keyword Arguments

- `output_file::String` or `output_files::Vector{String}: if specified, plot will
  be saved there. In methods 1 and 2, if not specified, output files will be
  computed using `get_forecast_filename`
- `title::String` or `titles::Vector{String}`
- `start_date::Nullable{Date}`
- `end_date::Nullable{Date}`
- `hist_label::String`
- `forecast_label::String`
- `hist_color::Colorant`
- `forecast_color::Colorant`
- `linestyle::Symbol`
- `bands_color::Colorant`
- `bands_pcts::Vector{String}`
- `bands_style::Symbol`: either `:fan` or `:line`
- `tick_size::Int`: x-axis (time) tick size in units of years
- `legend`
- `plot_handle::Plots.Plot`: a plot handle to add `history` and `forecast` to

**Method 1 only:**

- `forecast_string::String`
- `bdd_and_unbdd::Bool`: if true, then unbounded means and bounded bands are plotted

### Output

- `p::Plot`
"""
function plot_history_and_forecast(m::AbstractModel, var::Symbol, class::Symbol,
                                   input_type::Symbol, cond_type::Symbol;
                                   forecast_string::String = "",
                                   bdd_and_unbdd::Bool = false,
                                   fourquarter::Bool = false,
                                   output_file::String = "",
                                   title::String = "",
                                   kwargs...)

    plots = plot_history_and_forecast(m, [var], class, input_type, cond_type;
                                      forecast_string = forecast_string,
                                      bdd_and_unbdd = bdd_and_unbdd,
                                      fourquarter = fourquarter,
                                      output_files = isempty(output_file) ? String[] : [output_file],
                                      titles = isempty(title) ? String[] : [title],
                                      kwargs...)
    return plots[var]
end

function plot_history_and_forecast(m::AbstractModel, vars::Vector{Symbol}, class::Symbol,
                                   input_type::Symbol, cond_type::Symbol;
                                   forecast_string::String = "",
                                   bdd_and_unbdd::Bool = false,
                                   fourquarter::Bool = false,
                                   output_files::Vector{String} = String[],
                                   titles::Vector{String} = String[],
                                   kwargs...)
    # Read in MeansBands
    hist  = read_mb(m, input_type, cond_type, Symbol(fourquarter ? :hist4q : :hist, class),
                    forecast_string = forecast_string)
    fcast = read_mb(m, input_type, cond_type, Symbol(fourquarter ? :forecast4q : :forecast, class),
                    forecast_string = forecast_string, bdd_and_unbdd = bdd_and_unbdd)


    # Get output file names and titles if not provided
    if isempty(output_files)
        output_files = map(var -> get_forecast_filename(m, input_type, cond_type,
                                                        Symbol(fourquarter ? "forecast4q_" : "forecast_", var),
                                                        pathfcn = figurespath, fileformat = plot_extension()),
                           vars)
    end
    if isempty(titles)
        titles = map(var -> DSGE.describe_series(m, var, class), vars)
    end

    # Loop through variables
    plots = Dict{Symbol, Plots.Plot}()
    for (var, output_file, title) in zip(vars, output_files, titles)
        plots[var] = plot_history_and_forecast(var, hist, fcast;
                                               output_file = output_file, title = title,
                                               kwargs...)
    end

    return plots
end

function plot_history_and_forecast(var::Symbol, history::MeansBands, forecast::MeansBands;
                                   output_file::String = "",
                                   title::String = "",
                                   start_date::Nullable{Date} = Nullable{Date}(),
                                   end_date::Nullable{Date} = Nullable{Date}(),
                                   hist_label::String = "History",
                                   forecast_label::String = "Forecast",
                                   hist_color::Colorant = RGBA(0., 0., 0., 1.),
                                   forecast_color::Colorant = RGBA(1., 0., 0., 1.),
                                   linestyle::Symbol = :solid,
                                   bands_color::Colorant = RGBA(0., 0., 1., 0.1),
                                   bands_pcts::Vector{String} = String[],
                                   bands_style::Symbol = :fan,
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
    title!(p, title)

    # Plot bands
    if combined.metadata[:para] in [:full, :subset]
        bands_inds = get_bands_indices(var, history, forecast, hist_inds, fcast_inds)
        plot_bands!(p, var, combined, bands_style, bands_color,
                    linestyle = linestyle, pcts = bands_pcts, indices = bands_inds)
    end

    # Plot mean
    plot!(p, datenums[hist_inds],  combined.means[hist_inds,  var], label = hist_label,
          linewidth = 2, linestyle = linestyle, linecolor = hist_color)
    plot!(p, datenums[fcast_inds], combined.means[fcast_inds, var], label = forecast_label,
          linewidth = 2, linestyle = linestyle, linecolor = forecast_color)

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

"""
```
function plot_bands!(p, var, mb, style, color; linestyle = :solid,
    pcts = DSGE.which_density_bands(mb, uniquify = true), indices = Colon())
```

Plot `var` bands from `mb` on plot `p`. The `style` can be one of `:fan` or
`:line`, and the user can which bands to plot (`pcts`) and over which time
periods (`indices`).
"""
function plot_bands!(p::Plots.Plot, var::Symbol, mb::MeansBands,
                     style::Symbol, color::Colorant;
                     linestyle::Symbol = :solid,
                     pcts::Vector{String} = DSGE.which_density_bands(mb, uniquify = true),
                     indices = Colon())

    if isempty(pcts)
        pcts = which_density_bands(mb, uniquify = true)
    end

    datenums = map(quarter_date_to_number, mb.means[:date])

    if style == :fan
        for pct in pcts
            plot!(p, datenums[indices], mb.bands[var][indices, Symbol(pct, " UB")],
                  fillto = mb.bands[var][indices, Symbol(pct, " LB")],
                  label = "", color = color, α = 0.10)
        end

    elseif style == :line
        for pct in pcts
            plot!(p, datenums[indices], mb.bands[var][indices, Symbol(pct, " UB")],
                  label = "", color = color, linewidth = 2, linestyle = linestyle)
            plot!(p, datenums[indices], mb.bands[var][indices, Symbol(pct, " LB")],
                  label = "", color = color, linewidth = 2, linestyle = linestyle)
        end

    else
        error("Invalid style: " * string(style))
    end
end
