"""
```
plot_forecast_comparison(m_old, m_new, var, class, input_type, cond_type;
    title = "", kwargs...)

plot_forecast_comparison(m_old, m_new, vars, class, input_type, cond_type;
    forecast_string = "", bdd_and_unbdd = false, bands_pcts = [\"90.0%\"],
    old_hist_label = "", new_hist_label = "",
    old_forecast_label = \"Old Forecast\", new_forecast_label = \"New Forecast\",
    old_hist_color = :grey, new_hist_color = :black,
    old_forecast_color = :blue, new_forecast_color = :red,
    verbose = :low)
```

Plot forecasts from `m_old` and `m_new` of `var` or `vars`.

### Inputs

- `m_old::AbstractModel`
- `m_new::AbstractModel`
- `var::Symbol` or `vars::Vector{Symbol}`: e.g. `:obs_gdp` or `[:obs_gdp,
  :obs_nominalrate]`
- `class::Symbol`
- `input_type::Symbol`
- `cond_type::Symbol`

### Keyword Arguments

- `forecast_string::String`
- `bdd_and_unbdd::Bool`: if true, then unbounded means and bounded bands are plotted
- `bands_pcts::Vector{String}`: which bands to plot
- `old_hist_label::String`
- `new_hist_label::String`
- `old_forecast_label::String`
- `new_forecast_label::String`
- `old_hist_color::Colorant`
- `new_hist_color::Colorant`
- `old_forecast_color::Colorant`
- `new_forecast_color::Colorant`
- `plotroot::String`: if nonempty, plots will be saved in that directory
- `title::String` or `titles::Vector{String}`
- `verbose::Symbol`

### Output

- `p::Plot` or `plots::OrderedDict{Symbol, Plot}`
"""
function plot_forecast_comparison(m_old::AbstractModel, m_new::AbstractModel,
                                  var::Symbol, class::Symbol,
                                  input_type::Symbol, cond_type::Symbol;
                                  title::String = "",
                                  kwargs...)

    plots = plot_forecast_comparison(m_old, m_new, [var], class, input_type, cond_type;
                                     titles = isempty(title) ? String[] : [title],
                                     kwargs...)
    return plots[var]
end

function plot_forecast_comparison(m_old::AbstractModel, m_new::AbstractModel,
                                  vars::Vector{Symbol}, class::Symbol,
                                  input_type::Symbol, cond_type::Symbol;
                                  forecast_string::String = "",
                                  bdd_and_unbdd::Bool = false,
                                  old_hist_label::String = "",
                                  new_hist_label::String = "",
                                  old_forecast_label::String = "Old Forecast",
                                  new_forecast_label::String = "New Forecast",
                                  old_hist_color::Colorant = colorant"gray",
                                  new_hist_color::Colorant = colorant"black",
                                  old_forecast_color::Colorant = colorant"blue",
                                  new_forecast_color::Colorant = colorant"red",
                                  plotroot::String = "",
                                  titles::Vector{String} = String[],
                                  kwargs...)
    # Read in MeansBands
    histold  = read_mb(m_old, input_type, cond_type, Symbol(:hist, class),
                       forecast_string = forecast_string)
    forecastold = read_mb(m_old, input_type, cond_type, Symbol(:forecast, class),
                       forecast_string = forecast_string, bdd_and_unbdd = bdd_and_unbdd)
    histnew  = read_mb(m_new, input_type, cond_type, Symbol(:hist, class),
                       forecast_string = forecast_string)
    forecastnew = read_mb(m_new, input_type, cond_type, Symbol(:forecast, class),
                       forecast_string = forecast_string, bdd_and_unbdd = bdd_and_unbdd)

    # Get titles if not provided
    if isempty(titles)
        detexify_title = typeof(Plots.backend()) == Plots.GRBackend
        titles = map(var -> describe_series(m_new, var, class, detexify = detexify_title), vars)
    end

    # Set up common keyword arguments
    common_kwargs = Dict{Symbol, Any}()
    common_kwargs[:start_date]  = start_date
    common_kwargs[:end_date]    = end_date
    common_kwargs[:bands_pcts]  = bands_pcts
    common_kwargs[:bands_style] = :line

    # Loop through variables
    plots = OrderedDict{Symbol, Plots.Plot}()
    for (var, title) in zip(vars, titles)

        # Call recipe
        p = histforecast(var, histold, forecastold;
                         hist_label = old_hist_label, forecast_label = old_forecast_label,
                         hist_color = old_hist_color, forecast_color = old_forecast_color,
                         bands_color = old_forecast_color, linestyle = :solid,
                         common_kwargs, kwargs...)

        histforecast!(var, histnew, forecastnew;
                      hist_label = new_hist_label, forecast_label = new_forecast_label,
                      hist_color = new_hist_color, forecast_color = new_forecast_color,
                      bands_color = new_forecast_color, linestyle = :dash,
                      common_kwargs, kwargs...)

        # Save plot
        if !isempty(plotroot)
            output_file = joinpath(plotroot, "forecastcomp_" * detexify(string(var)) * "." *
                                   string(plot_extension()))
            save_plot(p, output_file, verbose = verbose)
        end
    end
    return plots
end
