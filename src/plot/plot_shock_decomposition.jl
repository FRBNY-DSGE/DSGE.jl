"""
```
plot_shock_decomposition(m, var, class, input_type, cond_type;
    forecast_string = "", groups = shock_groupings(m), output_file = "",
    title = "", kwargs...)

plot_shock_decomposition(m, vars, class, input_type, cond_type;
    forecast_string = "", groups = shock_groupings(m), output_files = [],
    titles = [], kwargs...)

plot_shock_decomposition(var, shockdec, trend, dettrend, hist, forecast, groups;
    output_file = "", title = "",
    hist_label = \"Detrended History\", forecast_label = \"Detrended Forecast\",
    hist_color = :black, forecast_color = :red, tick_size = 5, legend = :best)
```

Plot shock decomposition(s) for `var` or `vars`.

### Inputs

- `var::Symbol` or `vars::Vector{Symbol}`: variable(s) whose shock decomposition
  is to be plotted, e.g. `:obs_gdp` or `[:obs_gdp, :obs_nominalrate]`

**Methods 1 and 2 only:**

- `m::AbstractModel`
- `class::Symbol`
- `input_type::Symbol`
- `cond_type::Symbol`

**Method 3 only:**

- `shockdec::MeansBands`
- `trend::MeansBands`
- `dettrend::MeansBands`
- `hist::MeansBands`
- `forecast::MeansBands`
- `groups::Vector{ShockGroup}`

### Keyword Arguments

- `output_file::String` or `output_files::Vector{String}`: if specified, plot
  will be saved there. In methods 1 and 2, if not specified, output file(s) will
  be computed using `get_forecast_filename`
- `title::String` or `titles::Vector{String}`
- `hist_label::String`
- `forecast_label::String`
- `hist_color::Colorant`
- `forecast_color::Colorant`
- `tick_size::Int`: x-axis (time) tick size in units of years
- `legend`

**Methods 1 and 2 only:**

- `forecast_string::String`
- `groups::Vector{ShockGroup}`

### Output

- `p::Plot`
"""
function plot_shock_decomposition(m::AbstractModel, var::Symbol, class::Symbol,
                                  input_type::Symbol, cond_type::Symbol;
                                  forecast_string::String = "",
                                  groups::Vector{ShockGroup} = shock_groupings(m),
                                  output_file::String = "",
                                  title = "",
                                  kwargs...)

    plots = plot_shock_decomposition(m, [var], class, input_type, cond_type;
                                     forecast_string = forecast_string,
                                     groups = groups,
                                     output_files = isempty(output_file) ? String[] : [output_file],
                                     titles = isempty(title) ? String[] : [title],
                                     kwargs...)
    return plots[var]
end

function plot_shock_decomposition(m::AbstractModel, vars::Vector{Symbol}, class::Symbol,
                                  input_type::Symbol, cond_type::Symbol;
                                  forecast_string::String = "",
                                  groups::Vector{ShockGroup} = shock_groupings(m),
                                  output_files::Vector{String} = String[],
                                  titles::Vector{String} = String[],
                                  kwargs...)
    # Read in MeansBands
    output_vars = [Symbol(prod, class) for prod in [:shockdec, :trend, :dettrend, :hist, :forecast]]
    mbs = map(output_var -> read_mb(m, input_type, cond_type, output_var, forecast_string = forecast_string),
              output_vars)

    # Get output file names and titles if not provided
    if isempty(output_files)
        output_files = map(var -> get_forecast_filename(m, input_type, cond_type, Symbol("shockdec_", var),
                                                        pathfcn = figurespath, fileformat = plot_extension()),
                           vars)
    end
    if isempty(titles)
        titles = map(var -> describe_series(m, var, class), vars)
    end

    # Loop through variables
    plots = Dict{Symbol, Plots.Plot}()
    for (var, output_file, title) in zip(vars, output_files, titles)
        plots[var] = plot_shock_decomposition(var, mbs..., groups;
                                              output_file = output_file, title = title,
                                              ylabel = DSGE.series_ylabel(m, var, class),
                                              kwargs...)
    end
    return plots
end

function plot_shock_decomposition(var::Symbol, shockdec::MeansBands,
                                  trend::MeansBands, dettrend::MeansBands,
                                  hist::MeansBands, forecast::MeansBands,
                                  groups::Vector{ShockGroup};
                                  output_file::String = "",
                                  title = "",
                                  start_date::Date = shockdec.means[1, :date],
                                  end_date::Date = shockdec.means[end, :date],
                                  hist_label::String = "Detrended History",
                                  forecast_label::String = "Detrended Forecast",
                                  hist_color::Colorant = RGBA(0., 0., 0., 1.),
                                  forecast_color::Colorant = RGBA(1., 0., 0., 1.),
                                  tick_size::Int = 5,
                                  ylabel::String = "",
                                  legend = :best)

    # Construct DataFrame with detrended mean, deterministic trend, and all shocks
    df = prepare_means_table_shockdec(shockdec, trend, dettrend, var,
                                           mb_hist = hist, mb_forecast = forecast,
                                           detexify_shocks = false)

    # Dates
    start_ind, end_ind = get_date_limit_indices(start_date, end_date, df[:date])
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
    all_inds = start_ind:end_ind
    date_ticks = get_date_ticks(df[all_inds, :date], tick_size = tick_size)
    x0 = shockdec_date_to_x(quarter_number_to_date(date_ticks.start), df[start_ind, :date])
    x1 = shockdec_date_to_x(quarter_number_to_date(date_ticks.stop),  df[start_ind, :date])
    xstep = tick_size * 4
    x_ticks = x0:xstep:x1

    # Plot bars
    ngroups = length(groups)
    colors = map(x -> x.color, groups)
    labels = map(x -> x.name,  groups)
    cat_names = map(Symbol, labels)

    p = groupedbar(convert(Array, df[all_inds, cat_names]),
                   xtick = (x_ticks, date_ticks),
                   labels = reshape(labels, 1, ngroups),
                   color = reshape(colors, 1, ngroups),
                   linealpha = 0.0,
                   bar_width = 1.0,
                   legend = legend,
                   legendfont = Plots.Font("sans-serif", 5, :hcenter, :vcenter, 0.0, parse(Colorant, :black)),
                   title = title,
                   ylabel = ylabel)

    # Plot detrended mean
    hist_end_date = enddate_means(hist)
    hist_end_ind  = findfirst(df[:date], hist_end_date)

    plot!(p, df[start_ind:hist_end_ind, :x], df[start_ind:hist_end_ind, :detrendedMean],
          color = hist_color, linewidth = 2, label = hist_label, ylim = :auto)
    plot!(p, df[hist_end_ind:end_ind, :x], df[hist_end_ind:end_ind, :detrendedMean],
          color = forecast_color, linewidth = 2, label = forecast_label, ylim = :auto)

    # Save if output_file provided
    save_plot(p, output_file)

    return p
end
