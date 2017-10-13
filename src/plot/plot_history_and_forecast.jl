"""
```
plot_history_and_forecast(m, var, class, input_type, cond_type;
    title = "", kwargs...)

plot_history_and_forecast(m, vars, class, input_type, cond_type;
    forecast_string = "", bdd_and_unbdd = false,
    untrans = false, fourquarter = false,
    output_files = figurespath(m, \"forecast\"), titles = [],
    kwargs...)

plot_history_and_forecast(var, history, forecast; output_file = "",
    title = "", start_date = history.means[1, :date],
    end_date = forecast.means[end, :date],
    hist_label = \"History\", forecast_label = \"Forecast\",
    hist_color = :black, forecast_color = :red, linestyle = :solid,
    bands_color = :blue,
    bands_pcts = union(which_density_bands(history, uniquify = true),
                       which_density_bands(forecast, uniquify = true)),
    bands_style = :fan, label_bands = false, transparent_bands = true,
    tick_size = 5, ylabel = "", legend = :best, plot_handle = plot(),
    verbose = :low)
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
- `start_date::Date`
- `end_date::Date`
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
- `verbose::Symbol`

**Methods 1 and 2 only:**

- `forecast_string::String`
- `bdd_and_unbdd::Bool`: if true, then unbounded means and bounded bands are plotted
- `plotroot::String`: if nonempty, plots will be saved in that directory
- `untrans::Bool`: whether to plot untransformed (model units) history and forecast
- `fourquarter::Bool`: whether to plot four-quarter history and forecast

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
                                   untrans::Bool = false,
                                   fourquarter::Bool = false,
                                   plotroot::String = figurespath(m, "forecast"),
                                   titles::Vector{String} = String[],
                                   kwargs...)
    # Determine output_vars
    if untrans && fourquarter
        error("Only one of untrans or fourquarter can be true")
    elseif untrans
        hist_prod  = :histut
        fcast_prod = :forecastut
    elseif fourquarter
        hist_prod  = :hist4q
        fcast_prod = :forecast4q
    else
        hist_prod  = :hist
        fcast_prod = :forecast
    end

    # Read in MeansBands
    hist  = read_mb(m, input_type, cond_type, Symbol(hist_prod, class), forecast_string = forecast_string)
    fcast = read_mb(m, input_type, cond_type, Symbol(fcast_prod, class), forecast_string = forecast_string,
                    bdd_and_unbdd = bdd_and_unbdd)


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
                                  Symbol(fcast_prod, "_", detexify(var)),
                                  forecast_string = forecast_string,
                                  fileformat = plot_extension())
        end
        ylabel = series_ylabel(m, var, class, untrans = untrans, fourquarter = fourquarter)

        plots[var] = plot_history_and_forecast(var, hist, fcast;
                                               output_file = output_file, title = title,
                                               ylabel = ylabel,
                                               kwargs...)
    end
    return plots
end

function plot_history_and_forecast(var::Symbol, history::MeansBands, forecast::MeansBands;
                                   output_file::String = "",
                                   plot_handle::Plots.Plot = plot(),
                                   verbose::Symbol = :low,
                                   kwargs...)
    # Call recipe
    p = plot(plot_handle)
    histforecast!(var, history, forecast; kwargs...)

    # Save if output_file provided
    save_plot(p, output_file, verbose = verbose)

    return p
end

function plot_history_and_forecast(var::Symbol, histories::Vector{MeansBands}, forecasts::Vector{MeansBands};
                                   start_date::Date = histories[1].means[1, :date],
                                   end_date::Date = forecasts[1].means[end, :date],
                                   output_file::String = "",
                                   hist_label::Vector{String} = fill("History", length(histories)),
                                   forecast_label::Vector{String} = fill("Forecast", length(forecasts)),
                                   hist_color::Vector{Colorant} = Colorant[colorant"black" for i = 1:length(histories)],
                                   forecast_color::Vector{Colorant} = Colorant[colorant"red" for i = 1:length(forecasts)],
                                   bands_color::Vector{Colorant} = Colorant[colorant"blue" for i = 1:length(forecasts)],
                                   linestyle::Vector{Symbol} = fill(:solid, length(forecasts)),
                                   plot_handle::Plots.Plot = plot(),
                                   kwargs...)

    @assert length(histories) == length(forecasts) "histories and forecasts must be same length"

    for (i, (history, forecast)) in enumerate(zip(histories, forecasts))
        plot_handle = plot_history_and_forecast(var, history, forecast;
                            start_date = start_date, end_date = end_date,
                            output_file = output_file,
                            hist_label = hist_label[i], forecast_label = forecast_label[i],
                            hist_color = hist_color[i], forecast_color = forecast_color[i],
                            bands_color = bands_color[i],
                            linestyle = linestyle[i],
                            plot_handle = plot_handle,
                            kwargs...)
    end

    return plot_handle
end