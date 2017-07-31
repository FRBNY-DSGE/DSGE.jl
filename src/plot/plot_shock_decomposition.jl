"""
```
plot_shock_decomposition(m, var, class, input_type, cond_type; forecast_string = "",
    output_file = "", hist_label = \"Detrended History\", forecast_label = \"Detrended Forecast\",
    hist_color = :black, forecast_color = :red, tick_size = 5, legend = :best)

plot_shock_decomposition(var, shockdec, trend, dettrend, hist, forecast, groups;
    output_file = "", hist_label = \"Detrended History\", forecast_label = \"Detrended Forecast\",
    hist_color = :black, forecast_color = :red, tick_size = 5, legend = :best)
```

Plot shock decomposition for `var`.

### Inputs

- `var::Symbol`: variable whose shock decomposition is to be plotted, e.g. `:obs_gdp`

**Method 1 only:**

- `m::AbstractModel`
- `class::Symbol`
- `input_type::Symbol`
- `cond_type::Symbol`

**Method 2 only:**

- `shockdec::MeansBands`
- `trend::MeansBands`
- `dettrend::MeansBands`
- `hist::MeansBands`
- `forecast::MeansBands`
- `groups::Vector{ShockGroup}`

### Keyword Arguments

- `output_file::String`: if specified, plot will be saved there. In method 1, if
  `output_file` is not specified, one will be computed using `get_forecast_filename`
- `hist_label::String`
- `forecast_label::String`
- `hist_color::Colorant`
- `forecast_color::Colorant`
- `tick_size::Int`: x-axis (time) tick size in units of years
- `legend`

**Method 1 only:**

- `forecast_string::String`

### Output

- `p::Plot`
"""
function plot_shock_decomposition(m::AbstractModel, var::Symbol, class::Symbol,
                                  input_type::Symbol, cond_type::Symbol;
                                  forecast_string::String = "",
                                  output_file::String = "",
                                  hist_label::String = "Detrended History",
                                  forecast_label::String = "Detrended Forecast",
                                  hist_color::Colorant = RGBA(0., 0., 0., 1.),
                                  forecast_color::Colorant = RGBA(1., 0., 0., 1.),
                                  tick_size::Int = 5,
                                  legend = :best)
    # Read in MeansBands
    output_vars = [Symbol(prod, class) for prod in [:shockdec, :trend, :dettrend, :hist, :forecast]]
    mbs = map(output_var -> read_mb(m, input_type, cond_type, output_var, forecast_string = forecast_string),
              output_vars)

    # Get shock groupings
    groups = shock_groupings(m)

    # Get output file name if not provided
    if isempty(output_file)
        output_file = get_forecast_filename(m, input_type, cond_type, Symbol("shockdec_", var),
                                            pathfcn = figurespath, fileformat = plot_extension())
    end

    # Appeal to second method
    plot_shock_decomposition(var, mbs..., groups, output_file = output_file,
                             hist_label = hist_label, forecast_label = forecast_label,
                             hist_color = hist_color, forecast_color = forecast_color,
                             tick_size = tick_size, legend = legend)
end

function plot_shock_decomposition(var::Symbol, shockdec::MeansBands,
                                  trend::MeansBands, dettrend::MeansBands,
                                  hist::MeansBands, forecast::MeansBands,
                                  groups::Vector{ShockGroup};
                                  output_file::String = "",
                                  hist_label::String = "Detrended History",
                                  forecast_label::String = "Detrended Forecast",
                                  hist_color::Colorant = RGBA(0., 0., 0., 1.),
                                  forecast_color::Colorant = RGBA(1., 0., 0., 1.),
                                  tick_size::Int = 5,
                                  legend = :best)

    # Construct DataFrame with detrended mean, deterministic trend, and all shocks
    df = prepare_means_table_shockdec(shockdec, trend, dettrend, var,
                                           mb_hist = hist, mb_forecast = forecast,
                                           detexify_shocks = false)
    df[:datenum] = map(quarter_date_to_number, df[:date])
    df[:x] = map(date -> shockdec_date_to_x(date, df[1, :date]), df[:date])
    nperiods = size(df, 1)

    # Sum shock values for each group
    v0 = zeros(nperiods)
    for group in groups
        shock_vectors = [df[shock] for shock in group.shocks]
        shock_sum = reduce(+, v0, shock_vectors)
        df[Symbol(group.name)] = shock_sum
    end

    # x-axis ticks
    date_ticks = get_date_ticks(df[:date], tick_size = tick_size)
    x0 = shockdec_date_to_x(quarter_number_to_date(date_ticks.start), df[1, :date])
    x1 = shockdec_date_to_x(quarter_number_to_date(date_ticks.stop),  df[1, :date])
    xstep = tick_size * 4
    x_ticks = x0:xstep:x1

    # Plot bars
    ngroups = length(groups)
    colors = map(x -> x.color, groups)
    labels = map(x -> x.name,  groups)
    cat_names = map(Symbol, labels)

    p = groupedbar(convert(Array, df[cat_names]),
                   xtick = (x_ticks, date_ticks),
                   labels = reshape(labels, 1, ngroups),
                   color = reshape(colors, 1, ngroups),
                   linealpha = 0.0,
                   bar_width = 1.0,
                   legend = legend)

    # Plot detrended mean
    hist_end_date = enddate_means(hist)
    hist_end_ind  = findfirst(df[:date], hist_end_date)

    plot!(p, df[1:hist_end_ind, :x], df[1:hist_end_ind, :detrendedMean],
          color = hist_color, linewidth = 2, label = hist_label, ylim = :auto)
    plot!(p, df[hist_end_ind:end, :x], df[hist_end_ind:end, :detrendedMean],
          color = forecast_color, linewidth = 2, label = forecast_label, ylim = :auto)

    # Save if output_file provided
    save_plot(p, output_file)

    return p
end
