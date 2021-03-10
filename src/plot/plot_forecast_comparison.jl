"""
```
plot_forecast_comparison(m_old, m_new, var, class, input_type, cond_type;
    title = "", kwargs...)

plot_forecast_comparison(m_old, m_new, vars, class, input_type, cond_type;
    input_type_old = input_type, cond_type_old = cond_type,
    forecast_string = "", forecast_string_old = forecast_string,
    use_bdd = :unbdd, bands_pcts = [\"90.0%\"],
    old_names = Dict(:hist => "", :forecast => \"Old Forecast\"),
    new_names = Dict(:hist => "", :forecast => \"New Forecast\"),
    old_colors = Dict(:hist => :grey, :forecast => :blue, :bands => :blue),
    new_colors = Dict(:hist => :black, :forecast => :red, :bands => :red),
    old_alphas = Dict{Symbol, Float64}(),
    new_alphas = Dict{Symbol, Float64}(),
    old_styles = Dict(:hist => :dash, :forecast => :dash, :bands => :dash),
    new_styles = Dict(:hist => :solid, :forecast => :solid, :bands => :solid),
    plotroot = "", titles = [], verbose = :low)
```

Plot forecasts from `m_old` and `m_new` of `var` or `vars`.

### Inputs

- `m_old::AbstractDSGEModel`
- `m_new::AbstractDSGEModel`
- `var::Symbol` or `vars::Vector{Symbol}`: e.g. `:obs_gdp` or `[:obs_gdp,
  :obs_nominalrate]`
- `class::Symbol`
- `input_type::Symbol`
- `cond_type::Symbol`

### Keyword Arguments

- `input_type_old::Symbol`
- `cond_type_old::Symbol`
- `forecast_string::String`
- `forecast_string_old::String`
- `use_bdd::Symbol`: if :bdd_and_unbdd, then unbounded means and bounded bands are plotted
- `bands_pcts::Vector{String}`: which bands to plot
- `old_names::Dict{Symbol, String}`: maps keys `[:hist, :forecast, :bands]` to
  labels for old forecast
- `new_names::Dict{Symbol, String}`
- `old_colors::Dict{Symbol, Any}`: maps keys `[:hist, :forecast, :bands]` to
  colors for old forecast
- `new_colors::Dict{Symbol, Any}`
- `old_alphas::Dict{Symbol, Float64}`: maps keys `[:hist, :forecast, :bands]` to
  transparency values for old forecast
- `new_alphas::Dict{Symbol, Float64}`
- `old_styles::Dict{Symbol, Symbol}`: maps keys `[:hist, :forecast, :bands]` to
  linestyles for old forecast
- `new_styles::Dict{Symbol, Symbol}`
- `plotroot::String`: if nonempty, plots will be saved in that directory
- `title::String` or `titles::Vector{String}`
- `verbose::Symbol`

### Output

- `p::Plot` or `plots::OrderedDict{Symbol, Plot}`
"""
function plot_forecast_comparison(m_old::AbstractDSGEModel, m_new::AbstractDSGEModel,
                                  var::Symbol, class::Symbol,
                                  input_type::Symbol, cond_type::Symbol;
                                  title::String = "",
                                  kwargs...)

    plots = plot_forecast_comparison(m_old, m_new, [var], class, input_type, cond_type;
                                     titles = isempty(title) ? String[] : [title],
                                     kwargs...)
    return plots[var]
end

function plot_forecast_comparison(m_old::AbstractDSGEModel, m_new::AbstractDSGEModel,
                                  vars::Vector{Symbol}, class::Symbol,
                                  input_type::Symbol, cond_type::Symbol;
                                  input_type_old::Symbol = input_type,
                                  cond_type_old::Symbol = cond_type,
                                  forecast_string::String = "",
                                  forecast_string_old::String = forecast_string,
                                  use_bdd::Symbol = :unbdd,
                                  bands_pcts::Vector{String} = ["90.0%"],
                                  old_names = Dict(:hist => "", :forecast => "Old Forecast"),
                                  new_names = Dict(:hist => "", :forecast => "New Forecast"),
                                  old_colors = Dict(:hist => :grey, :forecast => :blue, :bands => :blue),
                                  new_colors = Dict(:hist => :black, :forecast => :red, :bands => :red),
                                  old_alphas = Dict{Symbol, Float64}(),
                                  new_alphas = Dict{Symbol, Float64}(),
                                  old_styles = Dict(:hist => :dash, :forecast => :dash, :bands => :dash),
                                  new_styles = Dict(:hist => :solid, :forecast => :solid, :bands => :solid),
                                  plotroot::String = "",
                                  titles::Vector{String} = String[],
                                  verbose::Symbol = :low,
                                  kwargs...)
    # Read in MeansBands
    histold = read_mb(m_old, input_type_old, cond_type_old, Symbol(:hist, class),
                      forecast_string = forecast_string_old)
    forecastold = read_mb(m_old, input_type_old, cond_type_old, Symbol(:forecast, class),
                       forecast_string = forecast_string_old, use_bdd = use_bdd)
    histnew = read_mb(m_new, input_type, cond_type, Symbol(:hist, class),
                      forecast_string = forecast_string)
    forecastnew = read_mb(m_new, input_type, cond_type, Symbol(:forecast, class),
                       forecast_string = forecast_string, use_bdd = use_bdd)

    # Get titles if not provided
    if isempty(titles)
        detexify_title = typeof(Plots.backend()) == Plots.GRBackend
        titles = map(var -> describe_series(m_new, var, class, detexify = detexify_title), vars)
    end

    # Loop through variables
    plots = OrderedDict{Symbol, Plots.Plot}()
    for (var, title) in zip(vars, titles)
        # Call recipe
        plots[var] = histforecast(var, histold, forecastold;
                                  names = old_names, colors = old_colors,
                                  alphas = old_alphas, styles = old_styles,
                                  bands_pcts = bands_pcts, bands_style = :line,
                                  title = title, ylabel = series_ylabel(m_new, var, class),
                                  kwargs...)

        histforecast!(var, histnew, forecastnew;
                      names = new_names, colors = new_colors,
                      alphas = new_alphas, styles = new_styles,
                      bands_pcts = bands_pcts, bands_style = :line, kwargs...)

        # Save plot
        if !isempty(plotroot)
            output_file = joinpath(plotroot, "forecastcomp_" * detexify(string(var)) * "." *
                                   string(plot_extension()))
            save_plot(plots[var], output_file, verbose = verbose)
        end
    end
    return plots
end
