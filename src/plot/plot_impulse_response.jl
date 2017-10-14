"""
```
plot_impulse_response(m, shock, var, class, input_type, cond_type;
    title = "", kwargs...)

plot_impulse_response(m, shock, vars, class, input_type, cond_type;
    forecast_string = "", plotroot = figurespath(m, \"forecast\"),
    titles = [], kwargs...)
```

Plot the responses of `var` to `shock`. By default, only 90% bands are plotted.

### Inputs

- `m::AbstractModel`
- `shock::Symbol`: e.g. `:g_sh`
- `var::Symbol` or `vars::Vector{Symbol}`: response variable(s), e.g. `:obs_gdp`
- `class::Symbol`
- `input_type::Symbol`
- `cond_type::Symbol`

### Keyword Arguments

- `forecast_string::String`
- `plotroot::String`: if nonempty, plots will be saved in that directory
- `title::String` or `titles::Vector{String}`
- `verbose::Symbol`

See `?irf` for additional keyword arguments, all of which can be passed
into `plot_history_and_forecast`.

### Output

- `p::Plot` or `plots::OrderedDict{Symbol, Plot}`
"""
function plot_impulse_response(m::AbstractModel, shock::Symbol, var::Symbol, class::Symbol,
                               input_type::Symbol, cond_type::Symbol;
                               title::String = "",
                               kwargs...)

    plots = plot_impulse_response(m, shock, [var], class, input_type, cond_type;
                                  titles = isempty(title) ? String[] : [title],
                                  kwargs...)
    return plots[var]
end

function plot_impulse_response(m::AbstractModel, shock::Symbol, vars::Vector{Symbol}, class::Symbol,
                               input_type::Symbol, cond_type::Symbol;
                               forecast_string::String = "",
                               plotroot::String = figurespath(m, "forecast"),
                               titles::Vector{String} = String[],
                               kwargs...)
    # Read in MeansBands
    mb = read_mb(m, input_type, cond_type, Symbol(:irf, class), forecast_string = forecast_string)

    # Get titles if not provided
    if isempty(titles)
        detexify_title = typeof(Plots.backend()) == Plots.GRBackend
        titles = map(var -> describe_series(m, var, class, detexify = detexify_title), vars)
    end

    # Loop through variables
    plots = OrderedDict{Symbol, Plots.Plot}()
    for (var, title) in zip(vars, titles)
        # Call recipe
        plots[var] = irf(shock, var, mb; title = title, kwargs...)

        # Save plot
        if !isempty(plotroot)
            output_file = get_forecast_filename(plotroot, filestring_base(m), input_type, cond_type,
                                                Symbol("irf_", shock, "_", detexify(var)),
                                                forecast_string = forecast_string,
                                                fileformat = plot_extension())
            save_plot(plots[var], output_file, verbose = verbose)
        end
    end
    return plots
end