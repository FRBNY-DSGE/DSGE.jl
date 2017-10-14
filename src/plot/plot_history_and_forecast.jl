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
    plot_handle = plot(), verbose = :low, kwargs...)
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

- `history::MeansBands`
- `forecast::MeansBands`

### Keyword Arguments

- `plot_handle::Plots.Plot`: a plot handle to add `history` and `forecast` to
- `title::String` or `titles::Vector{String}`
- `verbose::Symbol`

See `?histforecast` for additional keyword arguments, all of which can be passed
into `plot_history_and_forecast`.

**Methods 1 and 2 only:**

- `forecast_string::String`
- `bdd_and_unbdd::Bool`: if true, then unbounded means and bounded bands are plotted
- `plotroot::String`: if nonempty, plots will be saved in that directory
- `untrans::Bool`: whether to plot untransformed (model units) history and forecast
- `fourquarter::Bool`: whether to plot four-quarter history and forecast
- `title::String` or `titles::Vector{String}`

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