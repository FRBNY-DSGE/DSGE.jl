"""
```
plot_forecast_comparison(m_old, m_new, var, class, input_type, cond_type;
    forecast_string = "", bdd_and_unbdd = false, plotroot = "", title = "",
    kwargs...)

plot_forecast_comparison(m_old, m_new, vars, class, input_type, cond_type;
    forecast_string = "", bdd_and_unbdd = false, plotroot = "", titles = [],
    kwargs...)

plot_forecast_comparison(var, histold, fcastold, histnew, fcastnew;
    output_file = "", title = "", start_date = histold.means[1, :date],
    end_date = fcastnew.means[end, :date], bandpcts::Vector{String} = [\"90.0%\"],
    old_hist_label = "", new_hist_label = "",
    old_fcast_label = \"Old Forecast\", new_fcast_label = \"New Forecast\",
    old_hist_color = :grey, new_hist_color = :black,
    old_fcast_color = :blue, new_fcast_color = :red,
    tick_size = 2, ylabel = "", legend = :best)
```

### Inputs

- `var::Symbol` or `vars::Vector{Symbol}`: e.g. `:obs_gdp` or `[:obs_gdp,
  :obs_nominalrate]`

**Methods 1 and 2 only:**

- `m_old::AbstractModel`
- `m_new::AbstractModel`
- `class::Symbol`
- `input_type::Symbol`
- `output_type::Symbol`

**Method 3 only:**

- `histold::MeansBands`
- `fcastold::MeansBands`
- `histnew::MeansBands`
- `fcastnew::MeansBands`

### Keyword Arguments

- `title::String` or `titles::Vector{String}`
- `start_date::Date`
- `end_date::Date`
- `bandpcts::Vector{String}`: which bands to plot
- `old_hist_label::String`
- `new_hist_label::String`
- `old_fcast_label::String`
- `new_fcast_label::String`
- `old_hist_color::Colorant`
- `new_hist_color::Colorant`
- `old_fcast_color::Colorant`
- `new_fcast_color::Colorant`
- `tick_size::Int`: x-axis (time) tick size in units of years
- `ylabel::String`
- `legend`

**Methods 1 and 2 only:**

- `forecast_string::String`
- `bdd_and_unbdd::Bool`: if true, then unbounded means and bounded bands are plotted
- `plotroot::String`: if nonempty, plots will be saved in that directory

**Method 3 only:**

- `output_file::String`: if nonempty, plot will be saved in that path


### Output

- `p::Plot` or `plots::OrderedDict{Symbol, Plot}`
"""
function plot_forecast_comparison(m_old::AbstractModel, m_new::AbstractModel,
                                  var::Symbol, class::Symbol,
                                  input_type::Symbol, cond_type::Symbol;
                                  forecast_string::String = "",
                                  bdd_and_unbdd::Bool = false,
                                  plotroot::String = "",
                                  title::String = "",
                                  kwargs...)

    plots = plot_forecast_comparison(m_old, m_new, [var], class, input_type, cond_type;
                                     forecast_string = forecast_string,
                                     bdd_and_unbdd = bdd_and_unbdd,
                                     plotroot = plotroot,
                                     titles = isempty(title) ? String[] : [title],
                                     kwargs...)
    return plots[var]
end

function plot_forecast_comparison(m_old::AbstractModel, m_new::AbstractModel,
                                  vars::Vector{Symbol}, class::Symbol,
                                  input_type::Symbol, cond_type::Symbol;
                                  forecast_string::String = "",
                                  bdd_and_unbdd::Bool = false,
                                  plotroot::String = "",
                                  titles::Vector{String} = String[],
                                  kwargs...)
    # Read in MeansBands
    histold  = read_mb(m_old, input_type, cond_type, Symbol(:hist, class),
                       forecast_string = forecast_string)
    fcastold = read_mb(m_old, input_type, cond_type, Symbol(:forecast, class),
                       forecast_string = forecast_string, bdd_and_unbdd = bdd_and_unbdd)
    histnew  = read_mb(m_new, input_type, cond_type, Symbol(:hist, class),
                       forecast_string = forecast_string)
    fcastnew = read_mb(m_new, input_type, cond_type, Symbol(:forecast, class),
                       forecast_string = forecast_string, bdd_and_unbdd = bdd_and_unbdd)

    # Get titles if not provided
    if isempty(titles)
        detexify_title = typeof(Plots.backend()) == Plots.GRBackend
        titles = map(var -> DSGE.describe_series(m_new, var, class, detexify = detexify_title), vars)
    end

    # Loop through variables
    plots = OrderedDict{Symbol, Plots.Plot}()
    for (var, title) in zip(vars, titles)
        output_file = if isempty(plotroot)
            ""
        else
            joinpath(plotroot, "forecastcomp_" * string(var) * "." * plot_extension())
        end

        plots[var] = plot_forecast_comparison(var, histold, fcastold, histnew, fcastnew;
                                              output_file = output_file, title = title,
                                              ylabel = DSGE.series_ylabel(m_new, var, class),
                                              kwargs...)
    end
    return plots
end

function plot_forecast_comparison(var::Symbol,
                                  histold::MeansBands, fcastold::MeansBands,
                                  histnew::MeansBands, fcastnew::MeansBands;
                                  output_file::String = "",
                                  title::String = "",
                                  start_date::Date = histold.means[1, :date],
                                  end_date::Date = fcastnew.means[end, :date],
                                  bandpcts::Vector{String} = ["90.0%"],
                                  old_hist_label::String = "",
                                  new_hist_label::String = "",
                                  old_fcast_label::String = "Old Forecast",
                                  new_fcast_label::String = "New Forecast",
                                  old_hist_color::Colorant = parse(Colorant, :grey),
                                  new_hist_color::Colorant = parse(Colorant, :black),
                                  old_fcast_color::Colorant = parse(Colorant, :blue),
                                  new_fcast_color::Colorant = parse(Colorant, :red),
                                  tick_size::Int = 2,
                                  ylabel::String = "",
                                  legend = :best)

    # Initialize plot
    p = Plots.plot(legend = legend, title = title)

    # Set up common keyword arguments
    common_kwargs = Dict{Symbol, Any}()
    common_kwargs[:start_date]  = start_date
    common_kwargs[:end_date]    = end_date
    common_kwargs[:bands_pcts]  = bandpcts
    common_kwargs[:bands_style] = :line
    common_kwargs[:tick_size]   = tick_size
    common_kwargs[:legend]      = legend
    common_kwargs[:title]       = title
    common_kwargs[:ylabel]      = ylabel

    # Plot old and new histories/forecasts separately
    p = plot_history_and_forecast(var, histold, fcastold; plot_handle = p,
                                  hist_label = old_hist_label, forecast_label = old_fcast_label,
                                  hist_color = old_hist_color, forecast_color = old_fcast_color,
                                  bands_color = old_fcast_color, linestyle = :solid,
                                  common_kwargs...)
    p = plot_history_and_forecast(var, histnew, fcastnew; plot_handle = p,
                                  hist_label = new_hist_label, forecast_label = new_fcast_label,
                                  hist_color = new_hist_color, forecast_color = new_fcast_color,
                                  bands_color = new_fcast_color, linestyle = :dash,
                                  common_kwargs...)

    # Save if output_file provided
    save_plot(p, output_file)

    return p
end
