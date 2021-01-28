"""
```
plot_altpolicies(models, var, class, cond_type; kwargs...)

plot_altpolicies(models, vars, class, cond_type;
    forecast_string = "", altpol_string = "", untrans = false,
    fourquarter = false, plotroot = figurespath(m, \"forecast\"), titles = [],
    start_date = iterate_quarters(date_mainsample_end(models[1], -4)),
    end_date = iterate_quarters(date_forecast_start(models[1], 20)),
    kwargs...)
```

Plot `var` or `vars` forecasts under the alternative policies in `models`.

### Inputs

- `models::Vector{AbstractDSGEModel}`: vector of models, where each model has a
  different alternative policy
- `var::Symbol` or `vars::Vector{Symbol}`: variable(s) to be plotted,
  e.g. `:obs_gdp` or `[:obs_gdp, :obs_nominalrate]`
- `class::Symbol`
- `cond_type::Symbol`

### Keyword Arguments

- `forecast_string::String`
- `altpol_string::String`: identifies the group of alternative policies being
  plotted together. Required if `length(models) > 1`
- `start_date::Date`
- `end_date::Date`
- `untrans::Bool`: whether to plot untransformed forecast
- `fourquarter::Bool`: whether to plot four-quarter forecast
- `plotroot::String`: if nonempty, plots will be saved in that directory
- `title::String` or `titles::Vector{String}`

See `?histforecast` for additional keyword arguments, all of which can be passed
into `plot_altpolicies`.
"""
function plot_altpolicies(models::Vector{T}, var::Symbol, class::Symbol,
                          cond_type::Symbol; title::String = "",
                          kwargs...) where {T<:AbstractDSGEModel}
    plots = plot_altpolicies(models, [var], class, cond_type;
                             titles = isempty(title) ? String[] : [title],
                             kwargs...)
    return plots[var]
end

function plot_altpolicies(models::Vector{T}, vars::Vector{Symbol}, class::Symbol,
                          cond_type::Symbol;
                          forecast_string::String = "",
                          altpol_string::String = "",
                          start_date::Date = iterate_quarters(date_mainsample_end(models[1]), -4),
                          end_date::Date = iterate_quarters(date_forecast_start(models[1]), 20),
                          untrans::Bool = false,
                          fourquarter::Bool = false,
                          plotroot::String = figurespath(models[1], "forecast"),
                          titles::Vector{String} = String[],
                          verbose::Symbol = :low,
                          kwargs...) where {T<:AbstractDSGEModel}
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

    # Get titles if not provided
    if isempty(titles)
        detexify_title = typeof(Plots.backend()) == Plots.GRBackend
        titles = map(var -> describe_series(models[1], var, class, detexify = detexify_title), vars)
    end

    # Temporarily set models[1] altpolicy key to altpol_string for get_forecast_filename
    n_altpolicies = length(models)
    if n_altpolicies != 1
        if isempty(altpol_string)
            error("Must provide nonempty altpol_key if plotting multiple alternative policies")
        else
            altpol = alternative_policy(models[1])
            old_key = altpol.key
            altpol.key = Symbol(altpol_string)
        end
    end

    # Loop through variables
    plots = OrderedDict{Symbol, Plots.Plot}()
    for (var, title) in zip(vars, titles)
        # Loop through alternative policies
        plots[var] = plot()
        for (i, m) in enumerate(models)
            # Get AltPolicy
            altpolicy = alternative_policy(m)

            # Read in MeansBands
            hist     = read_mb(m, :mode, cond_type, Symbol(hist_prod, class), forecast_string = forecast_string)
            forecast = read_mb(m, :mode, cond_type, Symbol(fcast_prod, class), forecast_string = forecast_string)

            # Call recipe
            plots[var] = histforecast!(var, hist, forecast;
                                       start_date = start_date, end_date = end_date,
                                       names = Dict(:forecast =>string(altpolicy)),
                                       colors = Dict(:forecast => altpolicy.color),
                                       styles = Dict(:forecast => altpolicy.linestyle),
                                       ylabel = series_ylabel(m, var, class, untrans = untrans,
                                                              fourquarter = fourquarter),
                                       title = title, kwargs...)
        end

        # Save plot
        if !isempty(plotroot)
            output_file = get_forecast_filename(plotroot, filestring_base(models[1]), :mode, cond_type,
                                                Symbol("altpol", fcast_prod, "_", detexify(var)),
                                                forecast_string = forecast_string,
                                                fileformat = plot_extension())
            save_plot(plots[var], output_file, verbose = verbose)
        end
    end

    # Reset models[1] altpolicy key to original value
    if n_altpolicies != 1
        altpol.key = old_key
    end

    return plots
end
