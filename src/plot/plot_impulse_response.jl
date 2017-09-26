"""
```
plot_impulse_response(m, shock, var, class, input_type, cond_type;
    forecast_string = "", plotroot = figurespath(m, \"forecast\"),
    title = "", kwargs...)

plot_impulse_response(m, shock, vars, class, input_type, cond_type;
    forecast_string = "", plotroot = figurespath(m, \"forecast\"),
    titles = [], kwargs...)

plot_impulse_response(shock, var, mb; bands_pcts = [\"90.0%\"],
    output_file = "", title = "", flip = false, verbose = :low)
```

Plot the responses of `var` to `shock` from `mb`, possibly read in using
`read_mb` (depending on the method). By default, only 90% bands are plotted.

### Inputs

- `shock::Symbol`: e.g. `:g_sh`
- `var::Symbol` or `vars::Vector{Symbol}`: response variable(s), e.g. `:obs_gdp`

**Methods 1 and 2 only:**

- `m::AbstractModel`
- `class::Symbol`
- `input_type::Symbol`
- `cond_type::Symbol`

**Method 3 only:**

- `mb::MeansBands`

### Keyword Arguments

- `bands_pct::Vector{String}`
- `output_file::String`: if specified, plot will be saved there as a PDF
- `title::String` or `titles::Vector{String}`
- `flip::Bool`: whether to flip the sign of the impulse response while plotting
- `verbose::Symbol`

**Methods 1 and 2 only:**

- `forecast_string::String`
- `plotroot::String`: if nonempty, plots will be saved in that directory

**Method 3 only:**

- `output_file::String`: if nonempty, plot will be saved in that path

### Output

- `p::Plot` or `plots::OrderedDict{Symbol, Plot}`
"""
function plot_impulse_response(m::AbstractModel, shock::Symbol, var::Symbol, class::Symbol,
                               input_type::Symbol, cond_type::Symbol;
                               forecast_string::String = "",
                               plotroot::String = figurespath(m, "forecast"),
                               title::String = "",
                               kwargs...)

    plots = plot_impulse_response(m, shock, [var], class, input_type, cond_type;
                                  forecast_string = forecast_string,
                                  plotroot = plotroot,
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
        titles = map(var -> DSGE.describe_series(m, var, class, detexify = detexify_title), vars)
    end

    # Loop through variables
    plots = OrderedDict{Symbol, Plots.Plot}()
    for (var, title) in zip(vars, titles)
        output_file = if isempty(plotroot)
            ""
        else
            get_forecast_filename(plotroot, filestring_base(m), input_type, cond_type,
                                  Symbol("irf_", shock, "_", detexify(var)),
                                  forecast_string = forecast_string,
                                  fileformat = plot_extension())
        end

        plots[var] = plot_impulse_response(shock, var, mb;
                                           output_file = output_file, title = title,
                                           kwargs...)
    end
    return plots
end

function plot_impulse_response(shock::Symbol, var::Symbol, mb::MeansBands;
                               bands_pcts::Vector{String} = ["90.0%"],
                               output_file::String = "",
                               title::String = "",
                               flip::Bool = false,
                               verbose::Symbol = :low)

    # Plot title
    if isempty(title)
        title = string(var)
    end
    p = plot(title = title, margin = 10px)

    varshock = Symbol("$(var)__$(shock)")
    sign = flip ? -1 : 1

    # Plot bands
    for pct in bands_pcts
        plot!(p, sign * mb.bands[varshock][Symbol("$pct UB")],
              fillto = sign * mb.bands[varshock][Symbol("$pct LB")],
              label = "", color = :blue, α = 0.10)
    end

    # Plot mean
    plot!(p, sign * mb.means[varshock], label = "", linewidth = 2, linecolor = :black)

    # Save if `output_file` provided
    save_plot(p, output_file, verbose = verbose)

    return p
end
