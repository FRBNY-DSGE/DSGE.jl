"""
```
plot_forecast_comparison(m_old, m_new, var, class, input_type, cond_type;
    forecast_string = "", bdd_and_unbdd = false, output_file = "", title = "",
    kwargs...)

plot_forecast_comparison(m_old, m_new, vars, class, input_type, cond_type;
    forecast_string = "", bdd_and_unbdd = false, output_files = [], titles = [],
    kwargs...)

plot_forecast_comparison(var, histold, fcastold, histnew, fcastnew;
    output_file = "", title = "", start_date = Nullable{Date}(),
    end_date = Nullable{Date}(), bandpcts::Vector{String} = [\"90.0%\"],
    old_hist_label = "", new_hist_label = "",
    old_fcast_label = \"Old Forecast\", new_fcast_label = \"New Forecast\",
    old_hist_color = :grey, new_hist_color = :black,
    old_fcast_color = :blue, new_fcast_color = :red,
    tick_size = 2, legend = :best)
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

- `output_file::String` or `output_files::Vector{String}`: if specified, plot
  will be saved there as a PDF
- `title::String` or `titles::Vector{String}`
- `start_date::Nullable{Date}`
- `end_date::Nullable{Date}`
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
- `legend`

**Methods 1 and 2 only:**

- `forecast_string::String`
- `bdd_and_unbdd::Bool`: if true, then unbounded means and bounded bands are plotted

### Output

- `p::Plot`
"""
function plot_forecast_comparison(m_old::AbstractModel, m_new::AbstractModel,
                                  var::Symbol, class::Symbol,
                                  input_type::Symbol, cond_type::Symbol;
                                  forecast_string::String = "",
                                  bdd_and_unbdd::Bool = false,
                                  output_file::String = "",
                                  title::String = "",
                                  kwargs...)

    plot_forecast_comparison(m_old, m_new, [var], class, input_type, cond_type;
                             forecast_string = forecast_string,
                             bdd_and_unbdd = bdd_and_unbdd,
                             output_files = isempty(output_file) ? String[] : [output_file],
                             titles = isempty(title) ? String[] : [title],
                             kwargs...)
end

function plot_forecast_comparison(m_old::AbstractModel, m_new::AbstractModel,
                                  vars::Vector{Symbol}, class::Symbol,
                                  input_type::Symbol, cond_type::Symbol;
                                  forecast_string::String = "",
                                  bdd_and_unbdd::Bool = false,
                                  output_files::Vector{String} = String[],
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

    # Get output file names and titles if not provided
    if isempty(output_files)
        output_files = map(var -> "", vars)
    end
    if isempty(titles)
        titles = map(var -> describe_series(m_new, var, class), vars)
    end

    # Loop through variables
    for (var, output_file, title) in zip(vars, output_files, titles)

        plot_forecast_comparison(var, histold, fcastold, histnew, fcastnew;
                                 output_file = output_file, title = title,
                                 kwargs...)
    end
end

function plot_forecast_comparison(var::Symbol,
                                  histold::MeansBands, fcastold::MeansBands,
                                  histnew::MeansBands, fcastnew::MeansBands;
                                  output_file::String = "",
                                  title::String = "",
                                  start_date::Nullable{Date} = Nullable{Date}(),
                                  end_date::Nullable{Date} = Nullable{Date}(),
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
