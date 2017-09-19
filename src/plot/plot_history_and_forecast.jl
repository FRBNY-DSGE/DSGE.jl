"""
```
plot_history_and_forecast(m, var, class, input_type, cond_type;
    title = "", kwargs...)

plot_history_and_forecast(m, vars, class, input_type, cond_type;
    forecast_string = "", bdd_and_unbdd = false,
    output_files = figurespath(m, \"forecast\"), titles = [],
    kwargs...)

plot_history_and_forecast(var, history, forecast; output_file = "",
    title = "", start_date = Nullable{Date}(), end_date = Nullable{Date}(),
    hist_label = \"History\", forecast_label = \"Forecast\",
    hist_color = :black, forecast_color = :red, linestyle = :solid,
    bands_color = :blue,
    bands_pcts = union(which_density_bands(history, uniquify = true),
                       which_density_bands(forecast, uniquify = true)),
    bands_style = :fan, label_bands = false, transparent_bands = true,
    tick_size = 5, ylabel = "", legend = :best,
    plot_handle = plot())
```

Plot `var` or `vars` from `history` and `forecast`, possibly read in using
`read_mb` (depending on the method). If these correspond to a full-distribution
forecast, you can specify the `bands_style` and `bands_pcts`.

### Inputs

- `var::Symbol` or `vars::Vector{Symbol}`: variable(s) to be plotted,
  e.g. `:obs_gdp` or `[:obs_gdp, :obs_nominalrate]`

**Methods 1 and 2 only:**

- `m::AbstractModel`
- `class::Symbol`
- `input_type::Symbol`
- `cond_type::Symbol`

**Method 3 only:**

- `history::MeansBands` or `Vector{MeansBands}`
- `forecast::MeansBands` or `Vector{MeansBands}`

### Keyword Arguments

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
- `label_bands::Bool`
- `transparent_bands::Bool`
- `tick_size::Int`: x-axis (time) tick size in units of years
- `ylabel::String`
- `legend`
- `plot_handle::Plots.Plot`: a plot handle to add `history` and `forecast` to

**Methods 1 and 2 only:**

- `forecast_string::String`
- `bdd_and_unbdd::Bool`: if true, then unbounded means and bounded bands are plotted
- `plotroot::String`: if nonempty, plots will be saved in that directory

**Method 3 only:**

- `output_file::String`: if nonempty, plot will be saved in that path

### Output

- `p::Plot` or `plots::OrderedDict{Symbol, Plot}`
"""
function plot_history_and_forecast(m::AbstractModel, var::Symbol, class::Symbol,
                                   input_type::Symbol, cond_type::Symbol;
                                   title::String = "",
                                   kwargs...)

    plots = plot_history_and_forecast(m, [var], class, input_type, cond_type;
                                      titles = isempty(title) ? String[] : [title],
                                      kwargs...)
    return plots[var]
end

function plot_history_and_forecast(m::AbstractModel, vars::Vector{Symbol}, class::Symbol,
                                   input_type::Symbol, cond_type::Symbol;
                                   forecast_string::String = "",
                                   bdd_and_unbdd::Bool = false,
                                   fourquarter::Bool = false,
                                   plotroot::String = figurespath(m, "forecast"),
                                   titles::Vector{String} = String[],
                                   kwargs...)
    # Read in MeansBands
    hist  = read_mb(m, input_type, cond_type, Symbol(fourquarter ? :hist4q : :hist, class),
                    forecast_string = forecast_string)
    fcast = read_mb(m, input_type, cond_type, Symbol(fourquarter ? :forecast4q : :forecast, class),
                    forecast_string = forecast_string, bdd_and_unbdd = bdd_and_unbdd)


    # Get titles if not provided
    if isempty(titles)
        detexify_title = typeof(Plots.backend()) == Plots.GRBackend
        titles = map(var -> describe_series(m, var, class, detexify = detexify_title), vars)
    end

    # Loop through variables
    plots = OrderedDict{Symbol, Plots.Plot}()
    for (var, title) in zip(vars, titles)
        output_file = if isempty(plotroot)
            ""
        else
            get_forecast_filename(plotroot, filestring_base(m), input_type, cond_type,
                                  Symbol(fourquarter ? "forecast4q_" : "forecast_", detexify(var)),
                                  forecast_string = forecast_string,
                                  fileformat = plot_extension())
        end
        ylabel = series_ylabel(m, var, class, fourquarter = fourquarter)

        plots[var] = plot_history_and_forecast(var, hist, fcast;
                                               output_file = output_file, title = title,
                                               ylabel = ylabel,
                                               kwargs...)
    end
    return plots
end

function plot_history_and_forecast(var::Symbol, history::MeansBands, forecast::MeansBands;
                                   output_file::String = "",
                                   title::String = "",
                                   start_date::Date = history.means[1, :date],
                                   end_date::Date = forecast.means[end, :date],
                                   hist_label::String = "History",
                                   forecast_label::String = "Forecast",
                                   hist_color::Colorant = colorant"black",
                                   forecast_color::Colorant = colorant"red",
                                   linestyle::Symbol = :solid,
                                   bands_color::Colorant = colorant"blue",
                                   bands_pcts::Vector{String} = union(which_density_bands(history, uniquify = true),
                                                                      which_density_bands(forecast, uniquify = true)),
                                   bands_style::Symbol = :fan,
                                   label_bands::Bool = false,
                                   transparent_bands::Bool = true,
                                   tick_size::Int = 5,
                                   ylabel::String = "",
                                   legend = :best,
                                   plot_handle::Plots.Plot = plot(legend = legend))
    # Concatenate MeansBands
    combined = cat(history, forecast)

    # Dates
    start_ind, end_ind = get_date_limit_indices(start_date, end_date, combined.means[:date])
    datenums           = map(quarter_date_to_number, combined.means[:date])

    # Indices
    n_hist_periods  = size(history.means,  1)
    n_fcast_periods = size(forecast.means, 1)
    n_all_periods   = n_hist_periods + n_fcast_periods

    hist_inds  = start_ind:min(end_ind, n_hist_periods)
    fcast_inds = max(start_ind, n_hist_periods):end_ind
    all_inds   = start_ind:end_ind

    # Initialize plot
    p = plot_handle
    title!(p, title)
    yaxis!(p, ylabel = ylabel)

    # Plot bands
    if combined.metadata[:para] in [:full, :subset]
        bands_inds = get_bands_indices(var, history, forecast, hist_inds, fcast_inds)
        if !isempty(bands_inds)
            plot_bands!(p, var, combined, bands_style, bands_color,
                        linestyle = linestyle, pcts = bands_pcts, indices = bands_inds,
                        label_bands = label_bands, transparent_bands = transparent_bands)
        end
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
                                   plot_handle::Plots.Plot = plot(), kwargs...)

    @assert length(histories) == length(forecasts) "histories and forecasts must be same length"

    for (history, forecast) in zip(histories, forecasts)
       plot_handle = plot_history_and_forecast(var, history, forecast;
                                               plot_handle = plot_handle, kwargs...,)
    end

    return plot_handle
end

"""
```
function plot_bands!(p, var, mb, style, color; linestyle = :solid,
    pcts = which_density_bands(mb, uniquify = true), indices = Colon())
```

Plot `var` bands from `mb` on plot `p`. The `style` can be one of `:fan` or
`:line`, and the user can which bands to plot (`pcts`) and over which time
periods (`indices`).
"""
function plot_bands!(p::Plots.Plot, var::Symbol, mb::MeansBands,
                     style::Symbol, color::Colorant;
                     linestyle::Symbol = :solid,
                     pcts::Vector{String} = which_density_bands(mb, uniquify = true),
                     label_bands::Bool = false,
                     transparent_bands::Bool = true,
                     indices = Colon())

    datenums = map(quarter_date_to_number, mb.means[:date])

    # Sort percentages from largest to smallest
    sort!(pcts, rev = true)

    if style == :fan
        for (i, pct) in enumerate(pcts)
            # Determine label
            label = label_bands ? pct * " Bands" : ""

            # Determine color
            if transparent_bands
                pct_color = RGBA(color, 0.1*i)
                α = 0.1
            else
                pct_color = weighted_color_mean(0.1*i, color, colorant"white")
                α = 1.0
            end

            # Plot
            plot!(p, datenums[indices], mb.bands[var][indices, Symbol(pct, " UB")],
                  fillto = mb.bands[var][indices, Symbol(pct, " LB")],
                  label = label, color = pct_color, α = α)
        end

    elseif style == :line
        for pct in pcts
            # Determine labels
            upper_label = label_bands ? pct * " UB" : ""
            lower_label = label_bands ? pct * " LB" : ""

            # Plot
            plot!(p, datenums[indices], mb.bands[var][indices, Symbol(pct, " UB")],
                  label = upper_label, color = color, linewidth = 2, linestyle = linestyle)
            plot!(p, datenums[indices], mb.bands[var][indices, Symbol(pct, " LB")],
                  label = lower_label, color = color, linewidth = 2, linestyle = linestyle)
        end

    else
        error("Invalid style: " * string(style))
    end
end
