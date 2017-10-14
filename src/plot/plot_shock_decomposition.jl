"""
```
plot_shock_decomposition(m, var, class, input_type, cond_type;
    forecast_string = "", groups = shock_groupings(m),
    plotroot::String = figurespath(m, \"forecast\"), title = "",
    kwargs...)

plot_shock_decomposition(m, vars, class, input_type, cond_type;
    forecast_string = "", groups = shock_groupings(m),
    plotroot::String = figurespath(m, \"forecast\"), titles = [],
    kwargs...)

plot_shock_decomposition(var, shockdec, trend, dettrend, hist, forecast, groups;
    output_file = "", verbose = :low, kwargs...)
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

- `title::String` or `titles::Vector{String}`
- `verbose::Symbol`

See `?shockdec` for additional keyword arguments, all of which can be passed
into `plot_history_and_forecast`.

**Methods 1 and 2 only:**

- `forecast_string::String`
- `groups::Vector{ShockGroup}`
- `plotroot::String`: if nonempty, plots will be saved in that directory

**Method 3 only:**

- `output_file::String`: if nonempty, plot will be saved in that path

### Output

- `p::Plot` or `plots::OrderedDict{Symbol, Plot}`
"""
function plot_shock_decomposition(m::AbstractModel, var::Symbol, class::Symbol,
                                  input_type::Symbol, cond_type::Symbol;
                                  forecast_string::String = "",
                                  groups::Vector{ShockGroup} = shock_groupings(m),
                                  plotroot::String = figurespath(m, "forecast"),
                                  title = "",
                                  kwargs...)

    plots = plot_shock_decomposition(m, [var], class, input_type, cond_type;
                                     forecast_string = forecast_string,
                                     groups = groups,
                                     plotroot = plotroot,
                                     titles = isempty(title) ? String[] : [title],
                                     kwargs...)
    return plots[var]
end

function plot_shock_decomposition(m::AbstractModel, vars::Vector{Symbol}, class::Symbol,
                                  input_type::Symbol, cond_type::Symbol;
                                  forecast_string::String = "",
                                  groups::Vector{ShockGroup} = shock_groupings(m),
                                  plotroot::String = figurespath(m, "forecast"),
                                  titles::Vector{String} = String[],
                                  kwargs...)
    # Read in MeansBands
    output_vars = [Symbol(prod, class) for prod in [:shockdec, :trend, :dettrend, :hist, :forecast]]
    mbs = map(output_var -> read_mb(m, input_type, cond_type, output_var, forecast_string = forecast_string),
              output_vars)

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
                                  Symbol("shockdec_", detexify(var)),
                                  forecast_string = forecast_string,
                                  fileformat = plot_extension())
        end

        plots[var] = plot_shock_decomposition(var, mbs..., groups;
                                              output_file = output_file, title = title,
                                              ylabel = series_ylabel(m, var, class),
                                              kwargs...)
    end
    return plots
end

function plot_shock_decomposition(var::Symbol, shockdec::MeansBands,
                                  trend::MeansBands, dettrend::MeansBands,
                                  hist::MeansBands, forecast::MeansBands,
                                  groups::Vector{ShockGroup};
                                  output_file::String = "",
                                  verbose::Symbol = :low,
                                  kwargs...)

    # Call recipe
    shockdec(var, shockdec, trend, dettrend, hist, forecast, groups; kwargs...)

    # Save if output_file provided
    save_plot(p, output_file, verbose = verbose)

    return p
end