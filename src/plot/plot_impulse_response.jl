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
                               verbose::Symbol = :low,
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

@userplot Irf

"""
```
irf(shock, var, mb; flip = false, label_mean_bands = false,
    mean_color = :black, bands_color = :blue, bands_pcts = [\"90.0%\"])
```

User recipe called by `plot_impulse_response`.

### Inputs

- `shock::Symbol`: e.g. `:g_sh`
- `var::Symbol`: e.g. `:obs_gdp`
- `mb::MeansBands`

### Keyword Arguments

- `flip::Bool`: whether to flip the sign of the impulse response while plotting
- `label_mean_bands::Bool`
- `mean_color`
- `bands_color`
- `bands_pcts::Vector{String}`

Additionally, all Plots attributes (see docs.juliaplots.org/latest/attributes)
are supported as keyword arguments.
"""
irf

@recipe function f(irf::Irf;
                   flip = false,
                   label_mean_bands = false,
                   mean_color = :black,
                   bands_color = :blue,
                   bands_pcts = ["90.0%"])
    # Error checking
    if length(irf.args) != 3 || typeof(irf.args[1]) != Symbol || typeof(irf.args[2]) != Symbol ||
        typeof(irf.args[3]) != MeansBands

        error("irf must be given two Symbols and a MeansBands. Got $(typeof(irf.args))")
    end

    shock, var, mb = irf.args
    varshock = Symbol(var, "__", shock)
    sign = flip ? -1 : 1

    # Bands
    for pct in bands_pcts
        @series begin
            fillcolor := bands_color
            fillalpha := 0.1
            linealpha := 0
            label     := label_mean_bands ? "$pct Bands" : ""

            fillrange := sign * mb.bands[varshock][Symbol(pct, " UB")]
            sign * mb.bands[varshock][Symbol(pct, " LB")]
        end
    end

    # Mean
    @series begin
        label     := label_mean_bands ? "Mean" : ""
        linewidth := 2
        linecolor := mean_color
        sign * mb.means[varshock]
    end
end